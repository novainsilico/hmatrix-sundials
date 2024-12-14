{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE PatternSynonyms #-}
{-# LANGUAGE ViewPatterns #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

-- |
-- Solution of ordinary differential equation (ODE) initial value problems.
-- See <https://computation.llnl.gov/projects/sundials/sundials-software> for more detail.
module Numeric.Sundials.ARKode
  ( ARKMethod (..),
    solveC,
  )
where

import Control.Exception
import Control.Monad (when)
import Control.Monad.State
import Data.Maybe (fromMaybe)
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import Data.Void
import Foreign
import Foreign.C
import GHC.Generics
import GHC.Prim
import GHC.Stack
import Katip
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import Numeric.Sundials.MainLoop

-- | Available methods for ARKode
data ARKMethod
  = SDIRK_2_1_2
  | BILLINGTON_3_3_2
  | TRBDF2_3_3_2
  | KVAERNO_4_2_3
  | ARK324L2SA_DIRK_4_2_3
  | CASH_5_2_4
  | CASH_5_3_4
  | SDIRK_5_3_4
  | KVAERNO_5_3_4
  | ARK436L2SA_DIRK_6_3_4
  | KVAERNO_7_4_5
  | ARK548L2SA_DIRK_8_4_5
  | HEUN_EULER_2_1_2
  | BOGACKI_SHAMPINE_4_2_3
  | ARK324L2SA_ERK_4_2_3
  | ZONNEVELD_5_3_4
  | ARK436L2SA_ERK_6_3_4
  | SAYFY_ABURUB_6_3_4
  | CASH_KARP_6_4_5
  | FEHLBERG_6_4_5
  | DORMAND_PRINCE_7_4_5
  | ARK548L2SA_ERK_8_4_5
  | VERNER_8_5_6
  | FEHLBERG_13_7_8
  deriving (Eq, Ord, Show, Read, Generic, Bounded, Enum)

instance IsMethod ARKMethod where
  methodToInt SDIRK_2_1_2 = sDIRK_2_1_2
  methodToInt BILLINGTON_3_3_2 = bILLINGTON_3_3_2
  methodToInt TRBDF2_3_3_2 = tRBDF2_3_3_2
  methodToInt KVAERNO_4_2_3 = kVAERNO_4_2_3
  methodToInt ARK324L2SA_DIRK_4_2_3 = aRK324L2SA_DIRK_4_2_3
  methodToInt CASH_5_2_4 = cASH_5_2_4
  methodToInt CASH_5_3_4 = cASH_5_3_4
  methodToInt SDIRK_5_3_4 = sDIRK_5_3_4
  methodToInt KVAERNO_5_3_4 = kVAERNO_5_3_4
  methodToInt ARK436L2SA_DIRK_6_3_4 = aRK436L2SA_DIRK_6_3_4
  methodToInt KVAERNO_7_4_5 = kVAERNO_7_4_5
  methodToInt ARK548L2SA_DIRK_8_4_5 = aRK548L2SA_DIRK_8_4_5
  methodToInt HEUN_EULER_2_1_2 = hEUN_EULER_2_1_2
  methodToInt BOGACKI_SHAMPINE_4_2_3 = bOGACKI_SHAMPINE_4_2_3
  methodToInt ARK324L2SA_ERK_4_2_3 = aRK324L2SA_ERK_4_2_3
  methodToInt ZONNEVELD_5_3_4 = zONNEVELD_5_3_4
  methodToInt ARK436L2SA_ERK_6_3_4 = aRK436L2SA_ERK_6_3_4
  methodToInt SAYFY_ABURUB_6_3_4 = sAYFY_ABURUB_6_3_4
  methodToInt CASH_KARP_6_4_5 = cASH_KARP_6_4_5
  methodToInt FEHLBERG_6_4_5 = fEHLBERG_6_4_5
  methodToInt DORMAND_PRINCE_7_4_5 = dORMAND_PRINCE_7_4_5
  methodToInt ARK548L2SA_ERK_8_4_5 = aRK548L2SA_ERK_8_4_5
  methodToInt VERNER_8_5_6 = vERNER_8_5_6
  methodToInt FEHLBERG_13_7_8 = fEHLBERG_13_7_8

  methodType method =
    if methodToInt method < mIN_DIRK_NUM
      then Explicit
      else Implicit

foreign import ccall "wrapper"
  mkReport :: ReportErrorFnNew -> IO (FunPtr ReportErrorFnNew)

solveC :: Ptr CInt -> CConsts -> CVars (VS.MVector RealWorld) -> LogEnv -> IO CInt
solveC _ptrStop cconsts@CConsts {..} cvars@CVars {..} log_env =
  let report_error_new_api = wrapErrorNewApi (reportErrorWithKatip log_env)
      debug :: String -> StateT LoopState IO ()
      debug _s = do
        -- This SPAMs the logging system with a lot of informations, so we do
        -- not enable it by default
        -- Just comment in/out this line for now, we'll think about a nicer solution later
        -- liftIO (debugMsgWithKatip log_env s)
        pure ()
   in do
        -- Allocate the error reporting callback
        bracket (mkReport report_error_new_api) freeHaskellFunPtr $ \c_report_error -> do
          withSUNContext $ \sunctx -> do
            let init_loop =
                  ( LoopState
                      { -- /* input_ind tracks the current index into the c_sol_time array */
                        input_ind = 1,
                        -- /* output_ind tracks the current row into the c_output_mat matrix.
                        --    If differs from input_ind because of the extra rows corresponding to events. */
                        output_ind = 1,
                        -- /* event_ind tracks the current event number */
                        event_ind = 0,
                        -- /* t_start tracks the starting point of the integration in order to detect
                        --    empty integration interval and avoid a potential infinite loop;
                        --    see Note [CV_TOO_CLOSE]. Unlike T0, t_start is updated every time we
                        --    restart the solving after handling (or not) an event, or emitting
                        --    a requested time point.
                        --    Why not just look for the last recorded time in c_output_mat? Because
                        --    an event may have eventRecord = False and not be present there.
                        -- \*/
                        t_start = t0
                      }
                  )

                t0 = fromMaybe (error "no t0") $ c_sol_time VS.!? 0

            -- /* We need to update c_n_rows every time we update output_ind because
            --    of the possibility of early return (in which case we still need to assemble
            --    the partial results matrix). We could even work with c_n_rows only and ditch
            --    output_ind, but the inline-c expression is quite verbose, and output_ind is
            --    more convenient to use in index calculations.
            -- \*/
            VSM.write c_n_rows 0 (fromIntegral init_loop.output_ind)

            -- /* general problem parameters */

            -- /* Initialize data structures */

            -- /* Initialize odeMaxEventsReached to False */
            VSM.write c_diagnostics 10 0
            let implicit = c_method >= ARKODE_MIN_DIRK_NUM

            -- /* Create serial vector for solution */
            withNVector_Serial c_dim sunctx 6896 $ \y -> do
              -- /* Specify initial condition */
              VS.imapM_ (\i v -> cNV_Ith_S y i v) c_init_cond

              let withARKStep
                    | not implicit = \c_rhs -> withARKStepCreate c_rhs nullFunPtr
                    | otherwise = \c_rhs -> withARKStepCreate nullFunPtr c_rhs
              withARKStep c_rhs t0 y sunctx 8396 $ \cvode_mem -> do
                -- /* Set the error handler */
                cSUNContext_ClearErrHandlers sunctx >>= check 1093
                cSUNContext_PushErrHandler sunctx c_report_error nullPtr >>= check 1093

                when (c_fixedstep > 0.0) $ do
                  cARKodeSetFixedStep cvode_mem c_fixedstep

                -- /* Set the user data */
                cARKodeSetUserData cvode_mem c_rhs_userdata >>= check 1949

                -- /* Create serial vector for absolute tolerances */
                withNVector_Serial c_dim sunctx 6471 $ \tv -> do
                  -- /* Specify tolerances */
                  VS.imapM_ (\i v -> cNV_Ith_S tv i v) c_atol

                  cARKodeSetMinStep cvode_mem c_minstep >>= check 6433
                  cARKodeSetMaxNumSteps cvode_mem c_max_n_steps >>= check 9904
                  cARKodeSetMaxErrTestFails cvode_mem c_max_err_test_fails >>= check 2512

                  -- /* Specify the scalar relative tolerance and vector absolute tolerances */
                  cARKodeSVtolerances cvode_mem c_rtol tv >>= check 6212

                  -- /* Specify the root function */
                  cARKodeRootInit cvode_mem c_n_event_specs c_event_fn >>= check 6290
                  -- /* Disable the inactive roots warning; see https://git.novadiscovery.net/jinko/jinko/-/issues/2368 */
                  cARKodeSetNoInactiveRootWarn cvode_mem >>= check 6291

                  -- /* Initialize a jacobian matrix and solver */
                  let withLinearSolver f
                        | implicit = do
                            if (c_sparse_jac /= 0)
                              then do
                                withSUNSparseMatrix c_dim c_dim c_sparse_jac CSC_MAT sunctx 9061 $ \a -> do
                                  withSUNLinSol_KLU y a sunctx 9316 $ \ls -> do
                                    -- /* Attach matrix and linear solver */
                                    cARKodeSetLinearSolver cvode_mem ls a >>= check 2625
                                    f
                              else do
                                withSUNDenseMatrix c_dim c_dim sunctx 9316 $ \a -> do
                                  withSUNLinSol_Dense y a sunctx 9316 $ \ls -> do
                                    -- /* Attach matrix and linear solver */
                                    cARKodeSetLinearSolver cvode_mem ls a >>= check 2625
                                    f
                        | otherwise = f

                  withLinearSolver $ do
                    -- /* Set the initial step size if there is one */
                    when (c_init_step_size_set /= 0) $ do
                      --   /* FIXME: We could check if the initial step size is 0 */
                      --   /* or even NaN and then throw an error                 */
                      cARKodeSetInitStep cvode_mem c_init_step_size >>= check 4010

                    -- /* Set the Jacobian if there is one */
                    when (c_jac_set /= 0 && implicit) $ do
                      cARKodeSetJacFn cvode_mem c_jac >>= check 3124

                    -- /* Store initial conditions */
                    VSM.write c_output_mat (0 * (fromIntegral c_dim + 1) + 0) (c_sol_time VS.! 0)
                    let go j
                          | j == c_dim = pure ()
                          | otherwise = do
                              VSM.write c_output_mat (0 * fromIntegral (c_dim + 1) + (fromIntegral j + 1)) =<< cNV_Ith_S' y (fromIntegral j)
                              go (j + 1)
                    go 0

                    c_ontimepoint (fromIntegral init_loop.output_ind)

                    if implicit
                      then do
                        cARKStepSetTableNum cvode_mem c_method (-1) >>= check 26643
                      else do
                        cARKStepSetTableNum cvode_mem (-1) c_method >>= check 26643

                    let -- NOTE: this may be turned so the test is not done at each iteration
                        -- But we do not really care
                        reInit cvode_mem t y = do
                          if not implicit
                            then do
                              liftIO $ cARKStepReInit cvode_mem c_rhs nullFunPtr t y >>= check 1576
                            else do
                              liftIO $ cARKStepReInit cvode_mem nullFunPtr c_rhs t y >>= check 1576

                    let solver_interface =
                          SolverInterface
                            { ode = cARKodeEvolve,
                              odeGetEstLocalErrors = cARKodeGetEstLocalErrors,
                              odeGetErrWeights = cARKodeGetErrWeights,
                              odeGetRootInfo = cARKodeGetRootInfo,
                              odeReInit = reInit,
                              ode_NORMAL = ARK_NORMAL,
                              ode_ROOT_RETURN = ARK_ROOT_RETURN,
                              ode_TOO_CLOSE = ARK_TOO_CLOSE,
                              ode_SUCCESS = ARK_SUCCESS
                            }

                    resM <- try $ execStateT (mainLoop sunctx solver_interface cvode_mem cconsts cvars y debug) init_loop
                    case resM of
                      Left (ReturnCode c)
                        | c == fromIntegral ARK_SUCCESS -> pure ARK_SUCCESS
                        | otherwise -> pure $ (fromIntegral c)
                      Left (ReturnCodeWithMessage _message c)
                        | c == fromIntegral ARK_SUCCESS -> pure ARK_SUCCESS
                        | otherwise -> pure $ (fromIntegral c)
                      Right finalState -> end cvode_mem finalState
                      Left (Break finalState) -> end cvode_mem finalState
                      Left (Finish finalState) -> end cvode_mem finalState
  where
    end cvode_mem finalState = do
      -- /* The number of actual roots we found */
      VSM.write c_n_events 0 (fromIntegral finalState.event_ind)

      -- /* Get some final statistics on how the solve progressed */
      nst <- cvGet cARKodeGetNumSteps cvode_mem
      VSM.write c_diagnostics 0 (fromIntegral nst)

      nst_a <- cvGet cARKodeGetNumStepAttempts cvode_mem
      VSM.write c_diagnostics 1 (fromIntegral nst_a)

      (nfe, nfi) <- alloca $ \nfe -> do
        alloca $ \nfi -> do
          errCode <- cARKStepGetNumRhsEvals cvode_mem nfe nfi
          when (errCode /= ARK_SUCCESS) $ do
            error $ "Failure during nfe/nfi get"

          (,) <$> peek nfe <*> peek nfi

      VSM.write c_diagnostics 2 (fromIntegral nfe)
      VSM.write c_diagnostics 3 (fromIntegral nfi)

      nsetups <- cvGet cARKodeGetNumLinSolvSetups cvode_mem
      VSM.write c_diagnostics 4 (fromIntegral nsetups)

      netf <- cvGet cARKodeGetNumErrTestFails cvode_mem
      VSM.write c_diagnostics 5 (fromIntegral netf)

      nni <- cvGet cARKodeGetNumNonlinSolvIters cvode_mem
      VSM.write c_diagnostics 6 (fromIntegral nni)

      let implicit = c_method >= ARKODE_MIN_DIRK_NUM
      if implicit
        then do
          ncfn <- cvGet cARKodeGetNumNonlinSolvConvFails cvode_mem
          VSM.write c_diagnostics 7 (fromIntegral ncfn)

          nje <- cvGet cARKodeGetNumJacEvals cvode_mem
          VSM.write c_diagnostics 8 (fromIntegral nje)

          nfeLS <- cvGet cARKodeGetNumLinRhsEvals cvode_mem
          VSM.write c_diagnostics 9 (fromIntegral nfeLS)
        else do
          VSM.write c_diagnostics 7 0
          VSM.write c_diagnostics 8 0
          VSM.write c_diagnostics 9 0

      pure ARK_SUCCESS

--  |]

{- Note [CV_TOO_CLOSE]
   ~~~~~~~~~~~~~~~~~~~
   One edge condition that may occur is that an event time may exactly
   coincide with a solving time (e.g. they are both exactly equal to an
   integer). Then the following will happen:

   * Sundials will indicate a root at t1.
   * We will handle the event and re-initialize the system at t1.
   * We restart Sundials with the tout being equal to the next solving time,
     which also happens to be equal t1.
   * Sundials sees that the start and end solving times are equal, and
     returns the CV_TOO_CLOSE error.

   Calculating on our side when the start and end times are "too close" by
   Sundials standards is a bit complicated (see the code at the beginning
   of the cvHin function). It's much easier just to call Sundials and
   handle the error.

   For that, however, we need to make sure we ignore CV_TOO_CLOSE in our
   error handler so as not to confuse the end users with mysterious error
   messages in the logs.

   That said, we can't always rely on CV_TOO_CLOSE. When the initial step
   size is set, cvHin is not called, and CV_TOO_CLOSE is not triggered.
   Therefore we also add an explicit check to avoid an infinite loop of
   integrating over an empty interval.

   One exception to all of the above is a time-based event that may be
   scheduled for an exact same time as a time grid. In that case, we still
   handle it. We don't fall into an infinite loop because once we handle
   a time-based event, the next time-based event should be at a strictly
   later time.
-}

check :: (HasCallStack) => Int -> CInt -> IO ()
check retCode status
  | status == ARK_SUCCESS = pure ()
  | otherwise = throwIO (ReturnCode retCode)

-- | An opaque pointer to a ARKodeMem
newtype ARKodeMem = ARKodeMem (Ptr Void)
  deriving newtype (Storable)

withARKStepCreate ::
  (HasCallStack) =>
  ( FunPtr OdeRhsCType ->
    FunPtr OdeRhsCType ->
    CDouble ->
    N_Vector ->
    SUNContext ->
    Int ->
    (ARKodeMem -> IO c) ->
    IO c
  )
withARKStepCreate explicit implicit t0 y0 sunctx errCode f = do
  let create = do
        res@(ARKodeMem ptr) <- cARKStepCreate explicit implicit t0 y0 sunctx
        if ptr == nullPtr
          then throwIO $ ReturnCodeWithMessage "Error in cvodeCreate" errCode
          else pure res
      destroy p = do
        with p cARKStepFree
  bracket create destroy f

foreign import ccall "ARKStepCreate" cARKStepCreate :: FunPtr OdeRhsCType -> FunPtr OdeRhsCType -> CDouble -> N_Vector -> SUNContext -> IO ARKodeMem

foreign import ccall "ARKStepFree" cARKStepFree :: Ptr ARKodeMem -> IO ()

foreign import ccall "ARKodeSetUserData" cARKodeSetUserData :: ARKodeMem -> Ptr UserData -> IO CInt

foreign import ccall "ARKodeSetMinStep" cARKodeSetMinStep :: ARKodeMem -> CDouble -> IO CInt

foreign import ccall "ARKodeSetMaxNumSteps" cARKodeSetMaxNumSteps :: ARKodeMem -> SunIndexType -> IO CInt

foreign import ccall "ARKodeSetMaxErrTestFails" cARKodeSetMaxErrTestFails :: ARKodeMem -> CInt -> IO CInt

foreign import ccall "ARKodeSVtolerances" cARKodeSVtolerances :: ARKodeMem -> CDouble -> N_Vector -> IO CInt

foreign import ccall "ARKodeRootInit" cARKodeRootInit :: ARKodeMem -> CInt -> FunPtr EventConditionCType -> IO CInt

foreign import ccall "ARKodeSetNoInactiveRootWarn" cARKodeSetNoInactiveRootWarn :: ARKodeMem -> IO CInt

foreign import ccall "ARKodeSetLinearSolver" cARKodeSetLinearSolver :: ARKodeMem -> SUNLinearSolver -> SUNMatrix -> IO CInt

foreign import ccall "ARKodeSetInitStep" cARKodeSetInitStep :: ARKodeMem -> CDouble -> IO CInt

foreign import ccall "ARKodeEvolve" cARKodeEvolve :: ARKodeMem -> CDouble -> N_Vector -> Ptr CDouble -> CInt -> IO CInt

foreign import ccall "ARKStepReInit"
  cARKStepReInit ::
    ARKodeMem ->
    FunPtr OdeRhsCType ->
    FunPtr a2 ->
    CDouble ->
    N_Vector ->
    IO CInt

foreign import ccall "ARKodeGetRootInfo" cARKodeGetRootInfo :: ARKodeMem -> Ptr CInt -> IO CInt

foreign import ccall "ARKodeSetJacFn" cARKodeSetJacFn :: ARKodeMem -> FunPtr OdeJacobianCType -> IO CInt

foreign import ccall "ARKodeSetFixedStep" cARKodeSetFixedStep :: ARKodeMem -> CDouble -> IO ()

-- Note: the CInt are actually Enum and this could be enforced
foreign import ccall "ARKStepSetTableNum" cARKStepSetTableNum :: ARKodeMem -> CInt -> CInt -> IO CInt

foreign import ccall "ARKodeGetNumSteps" cARKodeGetNumSteps :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumLinSolvSetups" cARKodeGetNumLinSolvSetups :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumErrTestFails" cARKodeGetNumErrTestFails :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumNonlinSolvIters" cARKodeGetNumNonlinSolvIters :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumNonlinSolvConvFails" cARKodeGetNumNonlinSolvConvFails :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumJacEvals" cARKodeGetNumJacEvals :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKStepGetNumRhsEvals" cARKStepGetNumRhsEvals :: ARKodeMem -> Ptr CLong -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumLinRhsEvals" cARKodeGetNumLinRhsEvals :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumStepAttempts" cARKodeGetNumStepAttempts :: ARKodeMem -> Ptr CLong -> IO CInt

cvGet :: (HasCallStack) => (Storable b) => (ARKodeMem -> Ptr b -> IO CInt) -> ARKodeMem -> IO b
cvGet getter cvode_mem = do
  alloca $ \ptr -> do
    err <- getter cvode_mem ptr
    when (err /= ARK_SUCCESS) $ do
      error $ "Failure during cvGet"
    peek ptr

foreign import ccall "ARKodeGetEstLocalErrors" cARKodeGetEstLocalErrors :: ARKodeMem -> N_Vector -> IO CInt

foreign import ccall "ARKodeGetErrWeights" cARKodeGetErrWeights :: ARKodeMem -> N_Vector -> IO CInt
