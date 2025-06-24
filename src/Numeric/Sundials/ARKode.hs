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
import Data.Bool
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
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import Text.Printf (printf)

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

-- Tries to copy the previous semantic in C when we had:
data ReturnCode
  = -- return, with error code, skipping finish (diagnostics and cleanup)
    ReturnCode Int
  | -- return with log message, skipping finish (diags and cleanup)
    ReturnCodeWithMessage String Int
  | -- "goto" finish (do diag and cleanup)
    Finish LoopState
  | -- break the loop, hence go to finish (do diag and cleanup)
    Break LoopState
  deriving (Exception, Show)

-- | The loop state, most could be just carried by recursive function call, but
-- the C semantic was kinda interleaved, so doing that for now.
data LoopState = LoopState
  { -- output_ind tracks the current row into the c_output_mat matrix.
    -- if differs from input_ind because of the extra rows corresponding to events.
    output_ind :: Int,
    -- input_ind tracks the current index into the c_sol_time array
    input_ind :: Int,
    -- event_ind tracks the current event number
    event_ind :: Int,
    t_start :: CDouble
  }
  deriving (Show)

foreign import ccall "wrapper"
  mkReport :: ReportErrorFnNew -> IO (FunPtr ReportErrorFnNew)

solveC :: Ptr CInt -> CConsts -> CVars (VS.MVector RealWorld) -> LogEnv -> IO CInt
solveC _ptrStop CConsts {..} CVars {..} log_env =
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

              let withArkStep
                    | not implicit = \c_rhs -> withARKStepCreate c_rhs nullFunPtr
                    | otherwise = \c_rhs -> withARKStepCreate nullFunPtr c_rhs
              withArkStep c_rhs t0 y sunctx 8396 $ \cvode_mem -> do
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
                  case c_maxstep of
                    Just max_step -> cARKodeSetMaxStep cvode_mem max_step >>= check 6434
                    Nothing -> pure ()
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

                    let loop :: StateT LoopState IO ()
                        loop = do
                          s <- get
                          let ti = fromMaybe (error "Incorrect c_sol_time access") $ c_sol_time VS.!? s.input_ind

                          next_time_event <- liftIO c_next_time_event

                          --   // Haskell failure in the next time event function
                          when (next_time_event == -1) $ do
                            s <- get
                            liftIO $ throwIO $ Break s

                          when (next_time_event < s.t_start) $ do
                            s <- get
                            debug $ printf "time-based event is in the past: next event time = %.4f while we are at %.4f" (coerce next_time_event :: Double) (coerce s.t_start :: Double)
                            liftIO $ throwIO $ Finish s

                          let next_stop_time = min ti next_time_event
                          debug $ printf "Main loop iteration: t = %.17g (%a), next time point (ti) = %.17g, next time event = %.17g" (coerce s.t_start :: Double) (coerce ti :: Double) (coerce next_time_event :: Double)
                          (t, flag) <- liftIO $ alloca $ \t_ptr -> do
                            flag <- cARKodeEvolve cvode_mem next_stop_time y t_ptr ARK_NORMAL
                            t <- peek t_ptr
                            pure (t, flag)

                          debug $ printf "ARKode returned %d; now t = %.17g\n" (fromIntegral flag :: Int) (coerce t :: Double)
                          let root_based_event = flag == ARK_ROOT_RETURN
                          let time_based_event = t == next_time_event
                          (t, _flag) <-
                            if flag == ARK_TOO_CLOSE && not time_based_event
                              then do
                                --     /* See Note [CV_TOO_CLOSE]
                                --        No solving was required; just set the time t manually and continue
                                --        as if solving succeeded. */
                                debug $ printf "Got ARK_TOO_CLOSE; no solving was required; proceeding to t = %.17g" (coerce next_stop_time :: Double)
                                pure (next_stop_time, flag)
                              else do
                                s <- get
                                if t == next_stop_time && t == s.t_start && flag == ARK_ROOT_RETURN && not time_based_event
                                  then do
                                    --     /* See Note [CV_TOO_CLOSE]
                                    --        Probably the initial step size was set, and that's why we didn't
                                    --        get CV_TOO_CLOSE.
                                    --        Pretend that the root didn't happen, lest we keep handling it
                                    --        forever. */
                                    debug $ ("Got a root but t == t_start == next_stop_time; pretending it didn't happen" :: String)
                                    pure (t, ARK_SUCCESS)
                                  else do
                                    if not (flag == ARK_TOO_CLOSE && time_based_event) && flag < 0
                                      then do
                                        liftIO $ withNVector_Serial c_dim sunctx 12341234 $ \ele -> do
                                          liftIO $ withNVector_Serial c_dim sunctx 12341234 $ \weights -> do
                                            local_errors_flag <- liftIO $ cARKodeGetEstLocalErrors cvode_mem ele
                                            error_weights_flag <- liftIO $ cARKodeGetErrWeights cvode_mem weights
                                            when (local_errors_flag == ARK_SUCCESS && error_weights_flag == ARK_SUCCESS) $ do
                                              let go ix destination source
                                                    | ix == c_dim = pure ()
                                                    | otherwise = do
                                                        v <- peekElemOff source (fromIntegral ix)
                                                        VSM.write destination (fromIntegral ix) v
                                                        go (ix + 1) destination source
                                              go 0 c_local_error =<< cN_VGetArrayPointer ele
                                              go 0 c_var_weight =<< cN_VGetArrayPointer weights

                                              VSM.write c_local_error_set 0 1
                                            liftIO $ throwIO $ ReturnCode (fromIntegral flag)
                                      else pure (t, flag)

                          --   /* Store the results for Haskell */
                          s <- get
                          VSM.write c_output_mat (s.output_ind * (fromIntegral c_dim + 1) + 0) t
                          let go j
                                | j == c_dim = pure ()
                                | otherwise = do
                                    liftIO $ VSM.write c_output_mat (s.output_ind * fromIntegral (c_dim + 1) + (fromIntegral j + 1)) =<< cNV_Ith_S' y (fromIntegral j)
                                    go (j + 1)
                          go 0

                          s <- get
                          liftIO $ c_ontimepoint (fromIntegral s.output_ind)
                          modify $ \s -> s {output_ind = s.output_ind + 1}

                          s <- get
                          VSM.write c_n_rows 0 (fromIntegral s.output_ind)

                          when (root_based_event || time_based_event) $ do
                            debug ("Got an event")
                            when (fromIntegral s.event_ind >= c_max_events) $ do
                              debug ("Maximum number of events reached")
                              --       /* We reached the maximum number of events.
                              --          Either the maximum number of events is set to 0,
                              --          or there's a bug in our code below. In any case return an error.
                              --       */
                              liftIO $ throwIO (ReturnCode 8630)

                            --     /* How many events triggered? */
                            n_events_triggered <-
                              if not root_based_event
                                then pure 0
                                else do
                                  debug ("Handling root-based events")
                                  liftIO $ VSM.unsafeWith c_root_info $ \c_root_info_ptr -> do
                                    flag <- cARKodeGetRootInfo cvode_mem c_root_info_ptr
                                    when (flag < 0) $ do
                                      throwIO $ ReturnCode 2829
                                    let go i n_events_triggered
                                          | i >= c_n_event_specs = pure n_events_triggered
                                          | otherwise = do
                                              ev <- VSM.read c_root_info (fromIntegral i)
                                              let req_dir = c_requested_event_direction VS.! (fromIntegral i)

                                              if ev /= 0 && ev * req_dir >= 0
                                                then do
                                                  --           /* After the above call to ARKodeGetRootInfo, c_root_info has an
                                                  --           entry per EventSpec. Here we reuse the same array but convert it
                                                  --           into one that contains indices of triggered events. */
                                                  --           c_root_info[n_events_triggered++] = i;
                                                  VSM.write c_root_info n_events_triggered i
                                                  go (i + 1) (n_events_triggered + 1)
                                                else do
                                                  go (i + 1) n_events_triggered
                                    go 0 0

                            (record_events, stop_solver) <-
                              if (n_events_triggered > 0 || time_based_event)
                                then do
                                  debug $ printf "Calling the event handler; n_events_triggered = %d; time_based_event = %d" n_events_triggered (bool (0 :: Int) 1 time_based_event)
                                  (stop_solver, record_events, err) <- liftIO $ alloca $ \stop_solver_ptr -> alloca $ \record_event_ptr -> do
                                    --       /* Update the state with the supplied function */
                                    err <- VSM.unsafeWith c_root_info $ \c_root_info_ptr -> do
                                      err <- c_apply_event (fromIntegral n_events_triggered) c_root_info_ptr t (coerce y) (coerce y) stop_solver_ptr record_event_ptr
                                      pure err
                                    stop_solver <- peek stop_solver_ptr
                                    record_event <- peek record_event_ptr
                                    pure (stop_solver, record_event, err)

                                  --       // If the event handled failed internally, we stop the solving
                                  when (err /= 0) $ do
                                    s <- get
                                    liftIO $ throwIO $ Break s

                                  pure (record_events, stop_solver)
                                else pure (0, 0)

                            if record_events /= 0
                              then do
                                debug ("Recording events")
                                --       /* A corner case: if the time-based event triggers at the very beginning,
                                --          then we don't want to duplicate the initial row, so rewind it back.
                                --          Note that we do this only in the branch where record_events is true;
                                --          otherwise we may end up erasing the initial row (see below). */
                                s <- get
                                when (t == c_sol_time VS.! 0 && s.output_ind == 2) $ do
                                  --         /* c_n_rows will be updated below anyway */
                                  modify $ \s -> s {output_ind = s.output_ind - 1}

                                s <- get
                                VSM.write c_output_mat (s.output_ind * (fromIntegral c_dim + 1) + 0) t
                                let go j
                                      | j == c_dim = pure ()
                                      | otherwise = do
                                          liftIO $ VSM.write c_output_mat (s.output_ind * fromIntegral (c_dim + 1) + (fromIntegral j + 1)) =<< cNV_Ith_S' y (fromIntegral j)
                                          go (j + 1)
                                go 0

                                s <- get
                                liftIO $ c_ontimepoint $ fromIntegral s.output_ind
                                modify $ \s -> s {event_ind = s.event_ind + 1, output_ind = s.output_ind + 1}
                                s <- get
                                VSM.write c_n_rows 0 (fromIntegral s.output_ind)
                              else do
                                --       /* Remove the saved row — unless the event time also coincides with a requested time point */
                                when (t /= ti) $ do
                                  modify $ \s -> s {output_ind = s.output_ind - 1}
                                  s <- get
                                  VSM.write c_n_rows 0 (fromIntegral s.output_ind)
                            s <- get
                            stop_solver <-
                              if (fromIntegral s.event_ind >= c_max_events)
                                then do
                                  debug ("Reached max_events; returning")
                                  VSM.write c_diagnostics 10 1
                                  pure 1
                                else pure stop_solver
                            when (stop_solver /= 0) $ do
                              debug ("Stopping the hmatrix-sundials solver as requested")
                              s <- get
                              liftIO $ throwIO $ Finish s

                            when (n_events_triggered > 0 || time_based_event) $ do
                              debug ("Re-initializing the system")
                              if not implicit
                                then do
                                  liftIO $ cARKStepReInit cvode_mem c_rhs nullFunPtr t y >>= check 1576
                                else do
                                  liftIO $ cARKStepReInit cvode_mem nullFunPtr c_rhs t y >>= check 1576

                          when (t == ti) $ do
                            modify $ \s -> s {input_ind = s.input_ind + 1}
                            s <- get
                            when (s.input_ind >= fromIntegral c_n_sol_times) $ do
                              s <- get
                              liftIO $ throwIO $ Finish s

                          modify $ \s -> s {t_start = t}
                          loop
                    resM <- try $ execStateT loop init_loop
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

foreign import ccall "SUNContext_PushErrHandler" cSUNContext_PushErrHandler :: SUNContext -> FunPtr ReportErrorFnNew -> Ptr () -> IO CInt

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


{- NOTE [SAFETY]
   ~~~~~~~~~~~~~~~~~~~

All the call are marked (implicitly) as "safe" because they can, in theory, all call the log callback and hence "unsafe" is not relevant.

-}

-- An opaque pointer to a SUNContext
newtype SUNContext = SUNContext (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "SUNGetErrMsg" cSUNGetErrMsg :: CInt -> IO CString

withSUNContext :: (HasCallStack) => (SUNContext -> IO a) -> IO a
withSUNContext cont = do
  alloca $ \ptr -> do
    let create = do
          errCode <- cSUNContext_Create 0 ptr
          -- Exception safety here is a bit complex
          -- Let's assume that nothing will happen here
          when (errCode /= 0) $ do
            errMsg <- cSUNGetErrMsg errCode
            msg <- peekCString errMsg
            error msg

          -- Disable logging
          pure ()
        destroy = cSUNContext_Free ptr
    bracket_ create destroy $ do
      sunctx <- peek ptr

      -- We force disable the logging mecanism
      -- This stuff SPAMs stderr and generated exessive costs on our
      -- infrastructure
      alloca $ \loggerPtr -> do
        errCode <- cSUNContext_GetLogger sunctx loggerPtr
        when (loggerPtr == nullPtr || errCode /= 0) $ do
          error $ "logger is incomplete"
        logger <- peek loggerPtr

        withCString "/dev/null" $ \devNull -> do
          cSUNLogger_SetErrorFilename logger devNull
          cSUNLogger_SetWarningFilename logger devNull
          cSUNLogger_SetDebugFilename logger devNull
          cSUNLogger_SetInfoFilename logger devNull
          cont sunctx

check :: (HasCallStack) => Int -> CInt -> IO ()
check retCode status
  | status == ARK_SUCCESS = pure ()
  | otherwise = throwIO (ReturnCode retCode)

foreign import ccall "SUNContext_Create" cSUNContext_Create :: Int -> Ptr SUNContext -> IO CInt

foreign import ccall "SUNContext_Free" cSUNContext_Free :: Ptr SUNContext -> IO CInt

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

-- | An opaque pointer to an N_Vector
newtype N_Vector = N_Vector (Ptr Void)

foreign import ccall "N_VNew_Serial" unsafeN_VNew_Serial :: SunIndexType -> SUNContext -> IO N_Vector

foreign import ccall "N_VDestroy" unsafeN_VDestroy :: N_Vector -> IO ()

withNVector_Serial :: SunIndexType -> SUNContext -> Int -> (N_Vector -> IO a) -> IO a
withNVector_Serial t suncontext errCode f = do
  let create = do
        res@(N_Vector ptr) <- unsafeN_VNew_Serial t suncontext
        when (ptr == nullPtr) $ do
          throwIO $ ReturnCodeWithMessage "Failure in N_VNew_Serial" errCode
        pure res
  bracket create unsafeN_VDestroy f

cNV_Ith_S :: N_Vector -> Int -> CDouble -> IO ()
cNV_Ith_S (N_Vector ptr) i v = do
  qtr <- getContentPtr ptr
  rtr <- getData qtr
  pokeElemOff rtr i v

cNV_Ith_S' :: N_Vector -> Int -> IO CDouble
cNV_Ith_S' (N_Vector ptr) i = do
  qtr <- getContentPtr ptr
  rtr <- getData qtr
  peekElemOff rtr i

foreign import ccall "ARKodeSetUserData" cARKodeSetUserData :: ARKodeMem -> Ptr UserData -> IO CInt

foreign import ccall "ARKodeSetMinStep" cARKodeSetMinStep :: ARKodeMem -> CDouble -> IO CInt

foreign import ccall "ARKodeSetMaxStep" cARKodeSetMaxStep :: ARKodeMem -> CDouble -> IO CInt

foreign import ccall "ARKodeSetMaxNumSteps" cARKodeSetMaxNumSteps :: ARKodeMem -> SunIndexType -> IO CInt

foreign import ccall "ARKodeSetMaxErrTestFails" cARKodeSetMaxErrTestFails :: ARKodeMem -> CInt -> IO CInt

foreign import ccall "ARKodeSVtolerances" cARKodeSVtolerances :: ARKodeMem -> CDouble -> N_Vector -> IO CInt

foreign import ccall "ARKodeRootInit" cARKodeRootInit :: ARKodeMem -> CInt -> FunPtr EventConditionCType -> IO CInt

foreign import ccall "ARKodeSetNoInactiveRootWarn" cARKodeSetNoInactiveRootWarn :: ARKodeMem -> IO CInt

foreign import ccall "ARKodeSetLinearSolver" cARKodeSetLinearSolver :: ARKodeMem -> SUNLinearSolver -> SUNMatrix -> IO CInt

foreign import ccall "ARKodeSetInitStep" cARKodeSetInitStep :: ARKodeMem -> CDouble -> IO CInt

foreign import ccall "SUNSparseMatrix" cSUNSparseMatrix :: SunIndexType -> SunIndexType -> CInt -> CInt -> SUNContext -> IO SUNMatrix

foreign import ccall "SUNLinSol_KLU" cSUNLinSol_KLU :: N_Vector -> SUNMatrix -> SUNContext -> IO SUNLinearSolver

foreign import ccall "SUNDenseMatrix" cSUNDenseMatrix :: SunIndexType -> SunIndexType -> SUNContext -> IO SUNMatrix

foreign import ccall "SUNLinSol_Dense" cSUNLinSol_Dense :: N_Vector -> SUNMatrix -> SUNContext -> IO SUNLinearSolver

foreign import ccall "SUNMatDestroy" cSUNMatDestroy :: SUNMatrix -> IO ()

foreign import ccall "SUNLinSolFree" cSUNLinSolFree :: SUNLinearSolver -> IO CInt

withSUNDenseMatrix :: (HasCallStack) => SunIndexType -> SunIndexType -> SUNContext -> CInt -> (SUNMatrix -> IO a) -> IO a
withSUNDenseMatrix dim dim' sunctx errCode f = do
  let create = do
        mat@(SUNMatrix ptr) <- cSUNDenseMatrix dim dim' sunctx

        when (ptr == nullPtr) $ do
          throwIO $ ReturnCode $ fromIntegral errCode

        pure mat
  bracket create cSUNMatDestroy f

-- SUNLinSolFree ls

withSUNSparseMatrix :: (HasCallStack) => SunIndexType -> SunIndexType -> CInt -> CInt -> SUNContext -> CInt -> (SUNMatrix -> IO a) -> IO a
withSUNSparseMatrix dim dim' jac set sunctx errCode f = do
  let create = do
        mat@(SUNMatrix ptr) <- cSUNSparseMatrix dim dim' jac set sunctx

        when (ptr == nullPtr) $ do
          throwIO $ ReturnCode $ fromIntegral errCode

        pure mat
  bracket create cSUNMatDestroy f

withSUNLinSol_Dense :: (HasCallStack) => N_Vector -> SUNMatrix -> SUNContext -> CInt -> (SUNLinearSolver -> IO a) -> IO a
withSUNLinSol_Dense vec mat sunctx errCode f = do
  let create = do
        ls@(SUNLinearSolver ptr) <- cSUNLinSol_Dense vec mat sunctx

        when (ptr == nullPtr) $ do
          throwIO $ ReturnCode $ fromIntegral errCode

        pure ls
  bracket create cSUNLinSolFree f

withSUNLinSol_KLU :: (HasCallStack) => N_Vector -> SUNMatrix -> SUNContext -> CInt -> (SUNLinearSolver -> IO a) -> IO a
withSUNLinSol_KLU vec mat sunctx errCode f = do
  let create = do
        ls@(SUNLinearSolver ptr) <- cSUNLinSol_KLU vec mat sunctx

        when (ptr == nullPtr) $ do
          throwIO $ ReturnCode $ fromIntegral errCode

        pure ls
  bracket create cSUNLinSolFree f

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

-- | Opaque
newtype SUNMatrix = SUNMatrix (Ptr Void)
  deriving newtype (Storable)

newtype SUNLinearSolver = SUNLinearSolver (Ptr Void)
  deriving newtype (Storable)

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

foreign import ccall "SUNContext_ClearErrHandlers" cSUNContext_ClearErrHandlers :: SUNContext -> IO CInt

foreign import ccall "N_VGetArrayPointer" cN_VGetArrayPointer :: N_Vector -> IO (Ptr CDouble)

-- * Logs
newtype SUNLogger = SUNLogger (Ptr Void)
  deriving newtype Storable

foreign import ccall "SUNContext_GetLogger" cSUNContext_GetLogger :: SUNContext -> Ptr SUNLogger -> IO CInt

foreign import ccall "SUNLogger_SetErrorFilename" cSUNLogger_SetErrorFilename :: SUNLogger -> CString -> IO ()
foreign import ccall "SUNLogger_SetWarningFilename" cSUNLogger_SetWarningFilename :: SUNLogger -> CString -> IO ()
foreign import ccall "SUNLogger_SetInfoFilename" cSUNLogger_SetInfoFilename :: SUNLogger -> CString -> IO ()
foreign import ccall "SUNLogger_SetDebugFilename" cSUNLogger_SetDebugFilename :: SUNLogger -> CString -> IO ()


