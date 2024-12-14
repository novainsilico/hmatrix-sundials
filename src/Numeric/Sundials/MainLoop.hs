{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE OverloadedRecordDot #-}
{-# LANGUAGE PatternSynonyms #-}
{-# LANGUAGE ViewPatterns #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

module Numeric.Sundials.MainLoop where

import Control.Exception
import Control.Monad (when)
import Control.Monad.State
import Data.Bool
import Data.Maybe (fromMaybe)
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import Foreign
import Foreign.C
import GHC.Prim
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Text.Printf (printf)

data SolverInterface mem = SolverInterface
  {
    ode :: mem -> CDouble -> N_Vector -> Ptr CDouble -> CInt -> IO CInt,
    odeGetEstLocalErrors :: mem -> N_Vector -> IO CInt,
    odeGetErrWeights :: mem -> N_Vector -> IO CInt,
    odeGetRootInfo :: mem -> Ptr CInt -> IO CInt,
    odeReInit :: mem -> CDouble -> N_Vector -> IO (),
    ode_NORMAL :: CInt,
    ode_ROOT_RETURN :: CInt,
    ode_TOO_CLOSE :: CInt,
    ode_SUCCESS :: CInt
  }

mainLoop :: SUNContext -> SolverInterface mem -> mem -> CConsts -> CVars (VS.MVector RealWorld) -> N_Vector -> (String -> StateT LoopState IO ()) -> StateT LoopState IO ()
mainLoop sunctx SolverInterface {..} solver_mem CConsts {..} CVars {..} y debug = loop
  where
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
        flag <- ode solver_mem next_stop_time y t_ptr ode_NORMAL
        t <- peek t_ptr
        pure (t, flag)

      debug $ printf "CVode returned %d; now t = %.17g\n" (fromIntegral flag :: Int) (coerce t :: Double)
      let root_based_event = flag == ode_ROOT_RETURN
      let time_based_event = t == next_time_event
      (t, _flag) <-
        if flag == ode_TOO_CLOSE && not time_based_event
          then do
            --     /* See Note [ode_TOO_CLOSE]
            --        No solving was required; just set the time t manually and continue
            --        as if solving succeeded. */
            debug $ printf "Got ode_TOO_CLOSE; no solving was required; proceeding to t = %.17g" (coerce next_stop_time :: Double)
            pure (next_stop_time, flag)
          else do
            s <- get
            if t == next_stop_time && t == s.t_start && flag == ode_ROOT_RETURN && not time_based_event
              then do
                --     /* See Note [ode_TOO_CLOSE]
                --        Probably the initial step size was set, and that's why we didn't
                --        get ode_TOO_CLOSE.
                --        Pretend that the root didn't happen, lest we keep handling it
                --        forever. */
                debug $ ("Got a root but t == t_start == next_stop_time; pretending it didn't happen" :: String)
                pure (t, ode_SUCCESS)
              else do
                if not (flag == ode_TOO_CLOSE && time_based_event) && flag < 0
                  then do
                    liftIO $ withNVector_Serial c_dim sunctx 12341234 $ \ele -> do
                      liftIO $ withNVector_Serial c_dim sunctx 12341234 $ \weights -> do
                        flag <- liftIO $ odeGetEstLocalErrors solver_mem ele
                        flag' <- liftIO $ odeGetErrWeights solver_mem weights
                        when (flag == ode_SUCCESS && flag' == ode_SUCCESS) $ do
                          let go ix destination source
                                | ix == c_dim = pure ()
                                | otherwise = do
                                    v <- peekElemOff source (fromIntegral ix)
                                    VSM.write destination (fromIntegral ix) v
                                    go (ix + 1) destination source
                          go 0 c_local_error =<< cN_VGetArrayPointer ele
                          go 0 c_var_weight =<< cN_VGetArrayPointer weights

                          VSM.write c_local_error_set 0 1
                        liftIO $ throwIO (ReturnCode 45)
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
                flag <- odeGetRootInfo solver_mem c_root_info_ptr
                when (flag < 0) $ do
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
          liftIO $ odeReInit solver_mem t y

      when (t == ti) $ do
        modify $ \s -> s {input_ind = s.input_ind + 1}
        s <- get
        when (s.input_ind >= fromIntegral c_n_sol_times) $ do
          s <- get
          liftIO $ throwIO $ Finish s

      modify $ \s -> s {t_start = t}
      loop
