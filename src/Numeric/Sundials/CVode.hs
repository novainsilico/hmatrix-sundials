{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE PatternSynonyms #-}
{-# LANGUAGE ViewPatterns #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

-- | Solution of ordinary differential equation (ODE) initial value problems.
--
-- <https://computation.llnl.gov/projects/sundials/sundials-software>
module Numeric.Sundials.CVode
  ( CVMethod (..),
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

-- | Available methods for CVode
data CVMethod
  = ADAMS
  | BDF
  deriving (Eq, Ord, Show, Read, Generic, Bounded, Enum)

instance IsMethod CVMethod where
  methodToInt ADAMS = cV_ADAMS
  methodToInt BDF = cV_BDF
  methodType _ = Implicit

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

            -- /* Create serial vector for solution */
            withNVector_Serial c_dim sunctx 6896 $ \y -> do
              -- /* Specify initial condition */
              VS.imapM_ (\i v -> cNV_Ith_S y i v) c_init_cond

              -- // NB: Uses the Newton solver by default
              withCVodeMem c_method sunctx 8396 $ \cvode_mem -> do
                cCVodeInit cvode_mem c_rhs t0 y >>= check 1960

                -- /* Set the error handler */
                cSUNContext_ClearErrHandlers sunctx >>= check 1093
                cSUNContext_PushErrHandler sunctx c_report_error nullPtr >>= check 1093

                when (c_fixedstep > 0.0) $ do
                  throwIO $ ReturnCodeWithMessage "fixedStep cannot be used with CVode" 6426

                -- /* Set the user data */
                cCVodeSetUserData cvode_mem c_rhs_userdata >>= check 1949

                -- /* Create serial vector for absolute tolerances */
                withNVector_Serial c_dim sunctx 6471 $ \tv -> do
                  -- /* Specify tolerances */
                  VS.imapM_ (\i v -> cNV_Ith_S tv i v) c_atol

                  cCVodeSetMinStep cvode_mem c_minstep >>= check 6433
                  cCVodeSetMaxNumSteps cvode_mem c_max_n_steps >>= check 9904
                  cCVodeSetMaxErrTestFails cvode_mem c_max_err_test_fails >>= check 2512

                  -- /* Specify the scalar relative tolerance and vector absolute tolerances */
                  cCVodeSVtolerances cvode_mem c_rtol tv >>= check 6212

                  -- /* Specify the root function */
                  cCVodeRootInit cvode_mem c_n_event_specs c_event_fn >>= check 6290
                  -- /* Disable the inactive roots warning; see https://git.novadiscovery.net/jinko/jinko/-/issues/2368 */
                  cCVodeSetNoInactiveRootWarn cvode_mem >>= check 6291

                  -- /* Initialize a jacobian matrix and solver */
                  let withLinearSolver f = do
                        if (c_sparse_jac /= 0)
                          then do
                            withSUNSparseMatrix c_dim c_dim c_sparse_jac CSC_MAT sunctx 9061 $ \a -> do
                              withSUNLinSol_KLU y a sunctx 9316 $ \ls -> do
                                -- /* Attach matrix and linear solver */
                                cCVodeSetLinearSolver cvode_mem ls a >>= check 2625
                                f
                          else do
                            withSUNDenseMatrix c_dim c_dim sunctx 9316 $ \a -> do
                              withSUNLinSol_Dense y a sunctx 9316 $ \ls -> do
                                -- /* Attach matrix and linear solver */
                                cCVodeSetLinearSolver cvode_mem ls a >>= check 2625
                                f

                  withLinearSolver $ do
                    -- /* Set the initial step size if there is one */
                    when (c_init_step_size_set /= 0) $ do
                      --   /* FIXME: We could check if the initial step size is 0 */
                      --   /* or even NaN and then throw an error                 */
                      cCVodeSetInitStep cvode_mem c_init_step_size >>= check 4010

                    -- /* Set the Jacobian if there is one */
                    when (c_jac_set /= 0) $ do
                      cCVodeSetJacFn cvode_mem c_jac >>= check 3124

                    -- /* Store initial conditions */
                    VSM.write c_output_mat (0 * (fromIntegral c_dim + 1) + 0) (c_sol_time VS.! 0)
                    let go j
                          | j == c_dim = pure ()
                          | otherwise = do
                              VSM.write c_output_mat (0 * fromIntegral (c_dim + 1) + (fromIntegral j + 1)) =<< cNV_Ith_S' y (fromIntegral j)
                              go (j + 1)
                    go 0

                    c_ontimepoint (fromIntegral init_loop.output_ind)
                    let solver_interface =
                          SolverInterface
                            { ode = cCVode,
                              odeGetEstLocalErrors = cCVodeGetEstLocalErrors,
                              odeGetErrWeights = cCVodeGetErrWeights,
                              odeGetRootInfo = cCVodeGetRootInfo,
                              odeReInit = cCVodeReInit,
                              ode_NORMAL = CV_NORMAL,
                              ode_ROOT_RETURN = CV_ROOT_RETURN,
                              ode_TOO_CLOSE = CV_TOO_CLOSE,
                              ode_SUCCESS = CV_SUCCESS
                            }

                    resM <- try $ execStateT (mainLoop sunctx solver_interface cvode_mem cconsts cvars y debug) init_loop

                    case resM of
                      Left (ReturnCode c)
                        | c == fromIntegral CV_SUCCESS -> pure CV_SUCCESS
                        | otherwise -> pure $ (fromIntegral c)
                      Left (ReturnCodeWithMessage _message c)
                        | c == fromIntegral CV_SUCCESS -> pure CV_SUCCESS
                        | otherwise -> pure $ (fromIntegral c)
                      Right finalState -> end cvode_mem finalState
                      Left (Break finalState) -> end cvode_mem finalState
                      Left (Finish finalState) -> end cvode_mem finalState
  where
    end cvode_mem finalState = do
      -- /* The number of actual roots we found */
      VSM.write c_n_events 0 (fromIntegral finalState.event_ind)

      -- /* Get some final statistics on how the solve progressed */
      nst <- cvGet cCVodeGetNumSteps cvode_mem
      VSM.write c_diagnostics 0 (fromIntegral nst)

      -- /* FIXME */
      VSM.write c_diagnostics 1 0

      nfe <- cvGet cCVodeGetNumRhsEvals cvode_mem
      VSM.write c_diagnostics 2 (fromIntegral nfe)

      -- /* FIXME */
      VSM.write c_diagnostics 3 0

      nsetups <- cvGet cCVodeGetNumLinSolvSetups cvode_mem
      VSM.write c_diagnostics 4 (fromIntegral nsetups)

      netf <- cvGet cCVodeGetNumErrTestFails cvode_mem
      VSM.write c_diagnostics 5 (fromIntegral netf)

      nni <- cvGet cCVodeGetNumNonlinSolvIters cvode_mem
      VSM.write c_diagnostics 6 (fromIntegral nni)

      ncfn <- cvGet cCVodeGetNumNonlinSolvConvFails cvode_mem
      VSM.write c_diagnostics 7 (fromIntegral ncfn)

      nje <- cvGet cCVodeGetNumJacEvals cvode_mem
      VSM.write c_diagnostics 8 (fromIntegral nje)

      nfeLS <- cvGet cCVodeGetNumLinRhsEvals cvode_mem
      VSM.write c_diagnostics 9 (fromIntegral nfeLS)

      pure CV_SUCCESS

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
  | status == CV_SUCCESS = pure ()
  | otherwise = throwIO (ReturnCode retCode)

-- | An opaque pointer to a CVodeMem
newtype CVodeMem = CVodeMem (Ptr Void)
  deriving newtype (Storable)

withCVodeMem :: (HasCallStack) => CInt -> SUNContext -> Int -> (CVodeMem -> IO a) -> IO a
withCVodeMem method suncontext errCode f = do
  let create = do
        res@(CVodeMem ptr) <- cCVodeCreate method suncontext
        if ptr == nullPtr
          then throwIO $ ReturnCodeWithMessage "Error in cvodeCreate" errCode
          else pure res
      destroy p = do
        with p cCVodeFree
  bracket create destroy f

foreign import ccall "CVodeCreate" cCVodeCreate :: CInt -> SUNContext -> IO CVodeMem

foreign import ccall "CVodeFree" cCVodeFree :: Ptr CVodeMem -> IO ()

foreign import ccall "CVodeInit" cCVodeInit :: CVodeMem -> FunPtr OdeRhsCType -> SunRealType -> N_Vector -> IO CInt

foreign import ccall "CVodeSetUserData" cCVodeSetUserData :: CVodeMem -> Ptr UserData -> IO CInt

foreign import ccall "CVodeSetMinStep" cCVodeSetMinStep :: CVodeMem -> CDouble -> IO CInt

foreign import ccall "CVodeSetMaxNumSteps" cCVodeSetMaxNumSteps :: CVodeMem -> SunIndexType -> IO CInt

foreign import ccall "CVodeSetMaxErrTestFails" cCVodeSetMaxErrTestFails :: CVodeMem -> CInt -> IO CInt

foreign import ccall "CVodeSVtolerances" cCVodeSVtolerances :: CVodeMem -> CDouble -> N_Vector -> IO CInt

foreign import ccall "CVodeRootInit" cCVodeRootInit :: CVodeMem -> CInt -> FunPtr EventConditionCType -> IO CInt

foreign import ccall "CVodeSetNoInactiveRootWarn" cCVodeSetNoInactiveRootWarn :: CVodeMem -> IO CInt

foreign import ccall "CVodeSetLinearSolver" cCVodeSetLinearSolver :: CVodeMem -> SUNLinearSolver -> SUNMatrix -> IO CInt

foreign import ccall "CVodeSetInitStep" cCVodeSetInitStep :: CVodeMem -> CDouble -> IO CInt

foreign import ccall "CVode" cCVode :: CVodeMem -> CDouble -> N_Vector -> Ptr CDouble -> CInt -> IO CInt

foreign import ccall "CVodeReInit" cCVodeReInit :: CVodeMem -> CDouble -> N_Vector -> IO ()

foreign import ccall "CVodeGetRootInfo" cCVodeGetRootInfo :: CVodeMem -> Ptr CInt -> IO CInt

foreign import ccall "CVodeSetJacFn" cCVodeSetJacFn :: CVodeMem -> FunPtr OdeJacobianCType -> IO CInt

foreign import ccall "CVodeGetNumSteps" cCVodeGetNumSteps :: CVodeMem -> Ptr CLong -> IO CInt

foreign import ccall "CVodeGetNumLinSolvSetups" cCVodeGetNumLinSolvSetups :: CVodeMem -> Ptr CLong -> IO CInt

foreign import ccall "CVodeGetNumErrTestFails" cCVodeGetNumErrTestFails :: CVodeMem -> Ptr CLong -> IO CInt

foreign import ccall "CVodeGetNumNonlinSolvIters" cCVodeGetNumNonlinSolvIters :: CVodeMem -> Ptr CLong -> IO CInt

foreign import ccall "CVodeGetNumNonlinSolvConvFails" cCVodeGetNumNonlinSolvConvFails :: CVodeMem -> Ptr CLong -> IO CInt

foreign import ccall "CVodeGetNumJacEvals" cCVodeGetNumJacEvals :: CVodeMem -> Ptr CLong -> IO CInt

foreign import ccall "CVodeGetNumRhsEvals" cCVodeGetNumRhsEvals :: CVodeMem -> Ptr CLong -> IO CInt

foreign import ccall "CVodeGetNumLinRhsEvals" cCVodeGetNumLinRhsEvals :: CVodeMem -> Ptr CLong -> IO CInt

cvGet :: (HasCallStack) => (Storable b) => (CVodeMem -> Ptr b -> IO CInt) -> CVodeMem -> IO b
cvGet getter cvode_mem = do
  alloca $ \ptr -> do
    err <- getter cvode_mem ptr
    when (err /= CV_SUCCESS) $ do
      error $ "Failure during cvGet"
    peek ptr

foreign import ccall "CVodeGetEstLocalErrors" cCVodeGetEstLocalErrors :: CVodeMem -> N_Vector -> IO CInt

foreign import ccall "CVodeGetErrWeights" cCVodeGetErrWeights :: CVodeMem -> N_Vector -> IO CInt
