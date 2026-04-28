{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE PatternSynonyms #-}
{-# LANGUAGE ViewPatterns #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}
{-# LANGUAGE DataKinds #-}

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
import Data.Bool
import Data.Coerce (coerce)
import Data.Maybe (fromMaybe)
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import Foreign
import Foreign.C
import GHC.Generics
import GHC.Stack
import Katip
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import Text.Printf (printf)
import Numeric.Sundials.Bindings.CVode

-- | Available methods for CVode
data CVMethod
  = ADAMS
  | BDF
  deriving (Eq, Ord, Show, Read, Generic, Bounded, Enum)

instance IsMethod CVMethod where
  methodToInt ADAMS = cV_ADAMS
  methodToInt BDF = cV_BDF
  methodType _ = Implicit

solveC :: CConsts -> CVars (VS.MVector VSM.RealWorld) -> LogEnv -> IO (CInt, SundialsDiagnostics)
solveC CConsts {..} CVars {..} log_env =
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
                        --    see Note [#TOO_CLOSE]. Unlike T0, t_start is updated every time we
                        --    restart the solving after handling (or not) an event, or emitting
                        --    a requested time point.
                        --    Why not just look for the last recorded time in c_output_mat? Because
                        --    an event may have eventRecord = False and not be present there.
                        -- \*/
                        t_start = t0,
                        nb_reinit = 0,
                        max_events_reached = False,
                        current_diagnostics = mempty
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

            -- /* Create serial vector for solution */
            withNVector_Serial c_dim sunctx 6896 $ \y -> do
              -- /* Specify initial condition */
              VS.imapM_ (\i v -> cNV_Ith_S y i v) c_init_cond

              -- // NB: Uses the Newton solver by default
              withCVodeMem c_method sunctx 8396 $ \mem -> handleTermination @CVode #SUCCESS (getDiagnostics mem) $ do
                let getDiagnosticsCallback loopState = getDiagnostics mem loopState
                cCVodeInit mem c_rhs t0 y >>= check 1960

                -- /* Set the error handler */
                setErrorHandler sunctx c_report_error

                when (c_fixedstep > 0.0) $ do
                  throwIO $ ReturnCodeWithMessage "fixedStep cannot be used with CVode" 6426

                -- /* Set the user data */
                cCVodeSetUserData mem c_rhs_userdata >>= check 1949

                -- /* Create serial vector for absolute tolerances */
                withNVector_Serial c_dim sunctx 6471 $ \tv -> do
                  -- /* Specify tolerances */
                  VS.imapM_ (\i v -> cNV_Ith_S tv i v) c_atol

                  cCVodeSetMinStep mem c_minstep >>= check 6433
                  case c_maxstep of
                    Just max_step -> cCVodeSetMaxStep mem max_step >>= check 6434
                    Nothing -> pure ()
                  cCVodeSetMaxNumSteps mem c_max_n_steps >>= check 9904
                  cCVodeSetMaxErrTestFails mem c_max_err_test_fails >>= check 2512

                  -- /* Specify the scalar relative tolerance and vector absolute tolerances */
                  cCVodeSVtolerances mem c_rtol tv >>= check 6212

                  -- /* Specify the root function */
                  when (c_n_event_specs /= 0) $ do
                    cCVodeRootInit mem c_n_event_specs c_event_fn >>= check 6290

                    -- Set the root direction
                    VS.unsafeWith c_requested_event_direction $ \ptr -> do
                      cCVodeSetRootDirection mem ptr >>= check 5678909876
                    -- /* Disable the inactive roots warning; see https://git.novadiscovery.net/jinko/jinko/-/issues/2368 */
                    cCVodeSetNoInactiveRootWarn mem >>= check 6291

                  -- /* Initialize a jacobian matrix and solver */
                  let withLinearSolver f = do
                        if (c_sparse_jac /= 0)
                          then do
                            withSUNSparseMatrix c_dim c_dim c_sparse_jac CSC_MAT sunctx 9061 $ \a -> do
                              withSUNLinSol_KLU y a sunctx 9316 $ \ls -> do
                                -- /* Attach matrix and linear solver */
                                cCVodeSetLinearSolver mem ls a >>= check 2625
                                f
                          else do
                            withSUNDenseMatrix c_dim c_dim sunctx 9316 $ \a -> do
                              withSUNLinSol_Dense y a sunctx 9316 $ \ls -> do
                                -- /* Attach matrix and linear solver */
                                cCVodeSetLinearSolver mem ls a >>= check 2625
                                f

                  withLinearSolver $ do
                    -- /* Set the initial step size if there is one */
                    when (c_init_step_size_set /= 0) $ do
                      --   /* FIXME: We could check if the initial step size is 0 */
                      --   /* or even NaN and then throw an error                 */
                      cCVodeSetInitStep mem c_init_step_size >>= check 4010

                    -- /* Set the Jacobian if there is one */
                    when (c_jac_set /= 0) $ do
                      cCVodeSetJacFn mem c_jac >>= check 3124

                    -- /* Store initial conditions */
                    VSM.write c_output_mat (0 * (fromIntegral c_dim + 1) + 0) (c_sol_time VS.! 0)
                    let go j
                          | j == c_dim = pure ()
                          | otherwise = do
                              VSM.write c_output_mat (0 * fromIntegral (c_dim + 1) + (fromIntegral j + 1)) =<< cNV_Ith_S' y (fromIntegral j)
                              go (j + 1)
                    go 0

                    c_ontimepoint t0 c_output_mat 0 (getDiagnosticsCallback init_loop)

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
                            flag <- cCVode mem next_stop_time y t_ptr CV_NORMAL
                            -- When #TOO_CLOSE, the solver does not make any progress and does not update t_ptr
                            -- See https://github.com/llnl/sundials/issues/840
                            current_t <-
                              if flag == #TOO_CLOSE
                                then pure $ s.t_start
                                else peek t_ptr
                            pure (current_t, flag)

                          debug $ printf "CVode returned %d; now t = %.17g\n" (let Flag v = flag in (fromIntegral v :: Int)) (coerce t :: Double)
                          let root_based_event = flag == #ROOT_RETURN
                              time_based_event = t == next_time_event
                          t <-
                            if flag == #TOO_CLOSE && not time_based_event
                              then do
                                --     /* See Note [#TOO_CLOSE]
                                --        No solving was required; just set the time t manually and continue
                                --        as if solving succeeded. */
                                debug $ printf "Got #TOO_CLOSE; no solving was required; proceeding to t = %.17g" (coerce next_stop_time :: Double)
                                pure next_stop_time
                              else do
                                s <- get
                                if t == next_stop_time && t == s.t_start && flag == #ROOT_RETURN && not time_based_event
                                  then do
                                    --     /* See Note [#TOO_CLOSE]
                                    --        Probably the initial step size was set, and that's why we didn't
                                    --        get #TOO_CLOSE.
                                    --        Pretend that the root didn't happen, lest we keep handling it
                                    --        forever. */
                                    debug $ ("Got a root but t == t_start == next_stop_time; pretending it didn't happen" :: String)
                                    pure t
                                  else do
                                    if not (flag == #TOO_CLOSE && time_based_event) && isNegative flag
                                      then do
                                        liftIO $ withNVector_Serial c_dim sunctx 12341234 $ \ele -> do
                                          liftIO $ withNVector_Serial c_dim sunctx 12341234 $ \weights -> do
                                            local_errors_flag <- liftIO $ cCVodeGetEstLocalErrors mem ele
                                            error_weights_flag <- liftIO $ cCVodeGetErrWeights mem weights
                                            when (local_errors_flag == #SUCCESS && error_weights_flag == #SUCCESS) $ do
                                              let go ix destination source
                                                    | ix == c_dim = pure ()
                                                    | otherwise = do
                                                        v <- peekElemOff source (fromIntegral ix)
                                                        VSM.write destination (fromIntegral ix) v
                                                        go (ix + 1) destination source
                                              go 0 c_local_error =<< cN_VGetArrayPointer ele
                                              go 0 c_var_weight =<< cN_VGetArrayPointer weights

                                              VSM.write c_local_error_set 0 1
                                            liftIO $ throwIO $ ReturnCode (let Flag v = flag in fromIntegral v)
                                      else pure t

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
                          liftIO $ c_ontimepoint t c_output_mat (fromIntegral s.output_ind) (getDiagnosticsCallback s)
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
                                    flag <- cCVodeGetRootInfo mem c_root_info_ptr
                                    when (isNegative flag) $ do
                                      throwIO $ ReturnCode 2829
                                    let go i n_events_triggered
                                          | i >= c_n_event_specs = pure n_events_triggered
                                          | otherwise = do
                                              ev <- VSM.read c_root_info (fromIntegral i)
                                              let req_dir = c_requested_event_direction VS.! (fromIntegral i)

                                              if ev /= 0 && ev * req_dir >= 0
                                                then do
                                                  --           /* After the above call to CVodeGetRootInfo, c_root_info has an
                                                  --           entry per EventSpec. Here we reuse the same array but convert it
                                                  --           into one that contains indices of triggered events. */
                                                  --           c_root_info[n_events_triggered++] = i;
                                                  VSM.write c_root_info n_events_triggered i
                                                  go (i + 1) (n_events_triggered + 1)
                                                else do
                                                  go (i + 1) n_events_triggered
                                    go 0 0

                            (record_events, stop_solver, do_reinit) <-
                              if (n_events_triggered > 0 || time_based_event)
                                then do
                                  debug $ printf "Calling the event handler; n_events_triggered = %d; time_based_event = %d" n_events_triggered (bool (0 :: Int) 1 time_based_event)
                                  (stop_solver, record_events, err, do_reinit) <- liftIO $ alloca $ \stop_solver_ptr -> alloca $ \record_event_ptr -> alloca $ \do_reinit_ptr -> do
                                    --       /* Update the state with the supplied function */
                                    err <- VSM.unsafeWith c_root_info $ \c_root_info_ptr -> do
                                      err <- c_apply_event (fromIntegral n_events_triggered) c_root_info_ptr t (coerce y) (coerce y) nullPtr stop_solver_ptr record_event_ptr do_reinit_ptr (\x -> pure x)
                                      pure err
                                    stop_solver <- peek stop_solver_ptr
                                    record_event <- peek record_event_ptr
                                    do_reinit <- peek do_reinit_ptr
                                    pure (stop_solver, record_event, err, do_reinit)

                                  --       // If the event handled failed internally, we stop the solving
                                  when (err /= 0) $ do
                                    s <- get
                                    liftIO $ throwIO $ Break s

                                  pure (record_events, stop_solver, do_reinit)
                                else pure (0, 0, 0)

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
                                liftIO $ c_ontimepoint t c_output_mat (fromIntegral s.output_ind) (getDiagnosticsCallback s)
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
                                  modify $ \s -> s {max_events_reached = True}
                                  pure 1
                                else pure stop_solver
                            when (stop_solver /= 0) $ do
                              debug ("Stopping the hmatrix-sundials solver as requested")
                              s <- get
                              liftIO $ throwIO $ Finish s

                            when (do_reinit > 0) $ do
                              debug ("Re-initializing the system")

                              -- Before reinit, we get the diagnostics and update the current diagnostics
                              state <- get
                              diagnostics <- liftIO $ getDiagnostics mem state
                              put $ state {current_diagnostics = diagnostics}

                              -- This reset the diagnostics stored in cvode_mem
                              liftIO $ cCVodeReInit mem t y

                              modify $ \s -> s {nb_reinit = nb_reinit s + 1}

                          when (t == ti) $ do
                            modify $ \s -> s {input_ind = s.input_ind + 1}
                            s <- get
                            when (s.input_ind >= fromIntegral c_n_sol_times) $ do
                              s <- get
                              liftIO $ throwIO $ Finish s

                          modify $ \s -> s {t_start = t}
                          loop
                    execStateT loop init_loop

getDiagnostics :: CVodeMem -> LoopState -> IO SundialsDiagnostics
getDiagnostics cvode_mem loopState = do
  -- /* Get some final statistics on how the solve progressed */
  nst <- cvGet cCVodeGetNumSteps cvode_mem

  -- /* FIXME */
  let nst_a = 0 :: Int

  nfe <- cvGet cCVodeGetNumRhsEvals cvode_mem

  -- /* FIXME */
  let nfi = 0 :: Int

  nsetups <- cvGet cCVodeGetNumLinSolvSetups cvode_mem

  netf <- cvGet cCVodeGetNumErrTestFails cvode_mem

  nni <- cvGet cCVodeGetNumNonlinSolvIters cvode_mem

  ncfn <- cvGet cCVodeGetNumNonlinSolvConvFails cvode_mem

  nje <- cvGet cCVodeGetNumJacEvals cvode_mem

  nfeLS <- cvGet cCVodeGetNumLinRhsEvals cvode_mem

  gevals <- cvGet cCVodeGetNumGEvals cvode_mem

  let diagnostics =
        SundialsDiagnostics
          (fromIntegral $ nst)
          (fromIntegral $ nst_a)
          (fromIntegral $ nfe)
          (fromIntegral $ nfi)
          (fromIntegral $ nsetups)
          (fromIntegral $ netf)
          (fromIntegral $ nni)
          (fromIntegral $ ncfn)
          (fromIntegral $ nje)
          (fromIntegral $ nfeLS)
          False
          (fromIntegral gevals)
          0

  -- Some stats are not accumulated.
  pure $
    (diagnostics <> current_diagnostics loopState)
      { odeMaxEventsReached = loopState.max_events_reached,
        odeNumReinit = loopState.nb_reinit
      }

--  |]

{- Note [#TOO_CLOSE]
   ~~~~~~~~~~~~~~~~~~~
   One edge condition that may occur is that an event time may exactly
   coincide with a solving time (e.g. they are both exactly equal to an
   integer). Then the following will happen:

   * Sundials will indicate a root at t1.
   * We will handle the event and re-initialize the system at t1.
   * We restart Sundials with the tout being equal to the next solving time,
     which also happens to be equal t1.
   * Sundials sees that the start and end solving times are equal, and
     returns the #TOO_CLOSE error.

   Calculating on our side when the start and end times are "too close" by
   Sundials standards is a bit complicated (see the code at the beginning
   of the cvHin function). It's much easier just to call Sundials and
   handle the error.

   For that, however, we need to make sure we ignore #TOO_CLOSE in our
   error handler so as not to confuse the end users with mysterious error
   messages in the logs.

   That said, we can't always rely on #TOO_CLOSE. When the initial step
   size is set, cvHin is not called, and #TOO_CLOSE is not triggered.
   Therefore we also add an explicit check to avoid an infinite loop of
   integrating over an empty interval.

   One exception to all of the above is a time-based event that may be
   scheduled for an exact same time as a time grid. In that case, we still
   handle it. We don't fall into an infinite loop because once we handle
   a time-based event, the next time-based event should be at a strictly
   later time.
-}

check :: (HasCallStack) => Int -> Flag CVode -> IO ()
check retCode status
  | status == #SUCCESS = pure ()
  | otherwise = throwIO (ReturnCode retCode)

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

cvGet :: (HasCallStack) => (Storable b) => (CVodeMem -> Ptr b -> IO (Flag CVode)) -> CVodeMem -> IO b
cvGet getter cvode_mem = do
  alloca $ \ptr -> do
    err <- getter cvode_mem ptr
    when (err /= #SUCCESS) $ do
      error $ "Failure during cvGet"
    peek ptr
