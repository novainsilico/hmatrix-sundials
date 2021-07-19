-- | Common infrastructure for CVode/ARKode
{-# OPTIONS_GHC -Wno-name-shadowing #-}
module Numeric.Sundials.Common where

import Foreign.C.Types
import Foreign.Ptr
import Foreign.Storable (peek, poke)
import Foreign.C.String
import Numeric.Sundials.Foreign as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import Data.Maybe
import Data.Int
import Numeric.LinearAlgebra.HMatrix as H hiding (Vector)
import GHC.Prim
import Control.Monad.IO.Class
import Control.Monad.Cont
import Control.Exception
import Control.DeepSeq
import Katip
import Data.Aeson
import qualified Data.Text as T
import qualified Data.Text.Encoding as T
import qualified Data.ByteString as BS
import qualified Data.Map.Strict as Map
import Control.Monad.Reader
import GHC.Generics (Generic)
import Foreign.ForeignPtr
import Language.C.Types as CT
import Language.C.Inline.Context
import qualified Language.Haskell.TH as TH

-- | A collection of variables that we allocate on the Haskell side and
-- pass into the C code to be filled.
data CVars vec = CVars
  { c_diagnostics :: vec SunIndexType
    -- ^ Mutable vector to which we write diagnostic data while
    -- solving. Its size corresponds to the number of fields in
    -- 'SundialsDiagnostics'.
  , c_root_info :: vec CInt
    -- ^ Just a temporary vector (of the size equal to the number of event
    -- specs) that we use to get root info. Isn't used for output.
  , c_event_index :: vec CInt
    -- ^ For each event occurrence, this indicates which of the events
    -- occurred. Size: max_num_events.
  , c_event_time :: vec CDouble
    -- ^ For each event occurrence, this indicates the time of the
    -- occurrence. Size: max_num_events.
  , c_n_events :: vec CInt
    -- ^ Vector of size 1 that gives the total number of events occurred.
  , c_n_rows :: vec CInt
    -- ^ The total number of rows in the output matrix.
  , c_output_mat :: vec CDouble
    -- ^ The output matrix stored in the row-major order.
    -- Dimensions: (1 + dim) * (2 * max_events + nTs).
  , c_actual_event_direction :: vec CInt
    -- ^ Vector of size max_num_events that gives the direction of the
    -- occurred event.
  , c_local_error :: vec CDouble
    -- ^ Vector containing local error estimates. Size: the dimensionality
    -- of the system.
  , c_var_weight :: vec CDouble
    -- ^ Vector containing variable weights (derived from the tolerances).
    -- Size: the dimensionality of the system.
  , c_local_error_set :: vec CInt
    -- The flag (size 1) indicating whether c_local_error is filled with meaningful
    -- values. *Should be initialized with 0.*
  }

allocateCVars :: OdeProblem -> IO (CVars (VS.MVector RealWorld))
allocateCVars OdeProblem{..} = do
  let dim = VS.length odeInitCond
  c_diagnostics <- VSM.new 11
  c_root_info <- VSM.new $ V.length odeEventDirections
  c_event_index <- VSM.new odeMaxEvents
  c_event_time <- VSM.new odeMaxEvents
  c_actual_event_direction <- VSM.new odeMaxEvents
  c_n_events <- VSM.new 1
  c_n_rows <- VSM.new 1
  c_local_error <- VSM.new dim
  c_var_weight <- VSM.new dim
  c_local_error_set <- VSM.new 1
  VSM.write c_local_error_set 0 0
  c_output_mat <- VSM.new $
    (1 + dim) * (2 * odeMaxEvents + VS.length odeSolTimes)
  return CVars {..}

-- NB: the mutable CVars must not be used after this
freezeCVars :: CVars (VS.MVector RealWorld) -> IO (CVars VS.Vector)
freezeCVars CVars{..} = do
  c_diagnostics <- VS.unsafeFreeze c_diagnostics
  c_root_info <- VS.unsafeFreeze c_root_info
  c_event_index <- VS.unsafeFreeze c_event_index
  c_event_time <- VS.unsafeFreeze c_event_time
  c_actual_event_direction <- VS.unsafeFreeze c_actual_event_direction
  c_n_events <- VS.unsafeFreeze c_n_events
  c_n_rows <- VS.unsafeFreeze c_n_rows
  c_output_mat <- VS.unsafeFreeze c_output_mat
  c_local_error <- VS.unsafeFreeze c_local_error
  c_var_weight <- VS.unsafeFreeze c_var_weight
  c_local_error_set <- VS.unsafeFreeze c_local_error_set
  return CVars {..}

-- | Similar to 'CVars', except these are immutable values that are
-- accessed (read-only) by the C code and specify the system to be solved.
data CConsts = CConsts
  { c_dim :: SunIndexType -- ^ the dimensionality (number of variables/equations)
  , c_method :: CInt -- ^ the ODE method (specific to the solver)
  , c_n_sol_times :: CInt
  , c_sol_time :: VS.Vector CDouble
  , c_init_cond :: VS.Vector CDouble
  , c_rhs :: FunPtr OdeRhsCType
  , c_rhs_userdata :: Ptr UserData
  , c_rtol :: CDouble
  , c_atol :: VS.Vector CDouble
  , c_n_event_specs :: CInt
  , c_event_fn :: FunPtr EventConditionCType
  , c_apply_event
      :: CInt -- number of triggered events
      -> Ptr CInt -- event indices
      -> CDouble -- time
      -> Ptr T.SunVector -- y
      -> Ptr T.SunVector -- new y
      -> Ptr CInt -- (out) stop the solver?
      -> Ptr CInt -- (out) record the event?
      -> IO CInt
  , c_jac_set :: CInt
  , c_jac :: FunPtr OdeJacobianCType
  , c_sparse_jac :: CInt
      -- ^ If 0, use a dense matrix.
      -- If non-0, use a sparse matrix with that number of non-zero
      -- elements.
  , c_requested_event_direction :: VS.Vector CInt
  , c_next_time_event :: IO CDouble
  , c_max_events :: CInt
  , c_minstep :: CDouble
  , c_fixedstep :: CDouble
  , c_max_n_steps :: SunIndexType
  , c_max_err_test_fails :: CInt
  , c_init_step_size_set :: CInt
  , c_init_step_size :: CDouble
  }

data Solver = CVode | ARKode
  deriving Show

data MethodType = Explicit | Implicit
  deriving (Show, Eq)

class IsMethod method where
  methodToInt :: method -> CInt
  methodSolver :: Solver
  methodType :: method -> MethodType

withCConsts
  :: IsMethod method
  => ODEOpts method
  -> OdeProblem
  -> (CConsts -> IO r)
  -> IO r
withCConsts ODEOpts{..} OdeProblem{..} = runContT $ do
  let
    dim = VS.length c_init_cond
    c_init_cond = coerce odeInitCond
    c_dim = fromIntegral dim
    c_n_sol_times = fromIntegral . VS.length $ odeSolTimes
    c_sol_time = coerce odeSolTimes
    c_rtol = relTolerance odeTolerances
    c_atol = either (VS.replicate dim) id $ absTolerances odeTolerances
    c_minstep = coerce minStep
    c_fixedstep = coerce fixedStep
    c_max_n_steps = fromIntegral maxNumSteps
    c_max_err_test_fails = fromIntegral maxFail
    c_init_step_size_set = fromIntegral . fromEnum $ isJust initStep
    c_init_step_size = coerce . fromMaybe 0 $ initStep
    c_n_event_specs = fromIntegral $ V.length odeEventDirections
    c_requested_event_direction = V.convert $ V.map directionToInt odeEventDirections
    c_apply_event n_events event_indices_ptr t y_ptr y'_ptr stop_solver_ptr record_event_ptr = do
      event_indices <- vecFromPtr event_indices_ptr (fromIntegral n_events)
      -- Apparently there's no safe version of 
      y_vec <- peek y_ptr
      EventHandlerResult{..} <-
        odeEventHandler
          (coerce t :: Double)
          (coerce $ sunVecVals y_vec :: VS.Vector Double)
          (VS.map fromIntegral event_indices :: VS.Vector Int)
      poke y'_ptr $ SunVector
        { sunVecN = sunVecN y_vec
        , sunVecVals = coerce eventNewState
        }
      poke stop_solver_ptr . fromIntegral $ fromEnum eventStopSolver
      poke record_event_ptr . fromIntegral $ fromEnum eventRecord
      return 0
    c_max_events = fromIntegral odeMaxEvents
    c_next_time_event = coerce odeTimeBasedEvents
    c_jac_set = fromIntegral . fromEnum $ isJust odeJacobian
    c_sparse_jac = case jacobianRepr of
      SparseJacobian (T.SparsePattern spat) ->
        VS.sum (VS.map fromIntegral spat) +
        -- additionally, add diagonal zeros, as they'll be allocated too
        sum [ if spat VS.! (i + i * dim) == 0 then 1 else 0 | i <- [0 .. dim-1] ]
      DenseJacobian -> 0
    c_method = methodToInt odeMethod

  (c_rhs, c_rhs_userdata) <-
    case odeRhs of
      OdeRhsC ptr u -> return (ptr, u)
      OdeRhsHaskell fun -> do
        let
          funIO :: OdeRhsCType
          funIO t y f _ptr = do
            sv <- peek y
            r <- fun t (sunVecVals sv)
            poke f $ SunVector { sunVecN = sunVecN sv
                               , sunVecVals = r
                               }
            return 0
        funptr <- ContT $ bracket (mkOdeRhsC funIO) freeHaskellFunPtr
        return (funptr, nullPtr)
  c_jac <-
    case odeJacobian of
      Nothing   -> return nullFunPtr
      Just (OdeJacobianC fptr) -> return fptr
      Just (OdeJacobianHaskell jac_fn) -> do
      let
        funIO :: OdeJacobianCType
        funIO t y_ptr _fy_ptr jac_ptr _userdata _tmp1 _tmp2 _tmp3 = do
          y <- peek y_ptr
          let jac = matrixToSunMatrix $
                jac_fn
                  (coerce t :: Double)
                  (coerce $ sunVecVals y :: VS.Vector Double)
          case jacobianRepr of
            DenseJacobian -> poke jac_ptr jac
            SparseJacobian spat -> poke (castPtr jac_ptr) (T.SparseMatrix spat jac)
          return 0
      funptr <- ContT $ bracket (mkOdeJacobianC funIO) freeHaskellFunPtr
      return funptr
  c_event_fn <-
    case odeEventConditions of
      EventConditionsC fptr -> return fptr
      EventConditionsHaskell f -> do
      let
        funIO :: EventConditionCType
        funIO t y_ptr out_ptr _ptr = do
              y <- sunVecVals <$> peek y_ptr
              -- FIXME: We should be able to use poke somehow
              T.vectorToC (coerce f t y) (fromIntegral c_n_event_specs) out_ptr
              return 0
      funptr <- ContT $ bracket (mkEventConditionsC funIO) freeHaskellFunPtr
      return funptr
  return CConsts{..}

matrixToSunMatrix :: Matrix Double -> T.SunMatrix
matrixToSunMatrix m = T.SunMatrix { T.rows = nr, T.cols = nc, T.vals = vs }
  where
    nr = fromIntegral $ H.rows m
    nc = fromIntegral $ H.cols m
    vs = coerce . VS.concat $ toColumns m

-- Contrary to the documentation, it appears that CVodeGetRootInfo
-- may use both 1 and -1 to indicate a root, depending on the
-- direction of the sign change. See near the end of cvRootfind.
intToDirection :: Integral d => d -> Maybe CrossingDirection
intToDirection d =
  case d of
    1  -> Just Upwards
    -1 -> Just Downwards
    _  -> Nothing

-- | Almost inverse of 'intToDirection'. Map 'Upwards' to 1, 'Downwards' to
-- -1, and 'AnyDirection' to 0.
directionToInt :: Integral d => CrossingDirection -> d
directionToInt d =
  case d of
    Upwards -> 1
    Downwards -> -1
    AnyDirection -> 0

foreign import ccall "wrapper"
  mkOdeRhsC :: OdeRhsCType -> IO (FunPtr OdeRhsCType)

foreign import ccall "wrapper"
  mkOdeJacobianC :: OdeJacobianCType -> IO (FunPtr OdeJacobianCType)

foreign import ccall "wrapper"
  mkEventConditionsC :: EventConditionCType -> IO (FunPtr EventConditionCType)

assembleSolverResult
  :: OdeProblem
  -> CInt
  -> CVars VS.Vector
  -> IO (Either ErrorDiagnostics SundialsSolution)
assembleSolverResult OdeProblem{..} ret CVars{..} = do
  let
    dim = VS.length odeInitCond
    n_rows = fromIntegral . VS.head $ c_n_rows
    output_mat = coerce . reshape (dim + 1) . subVector 0 ((dim + 1) * n_rows) $ c_output_mat
    (local_errors, var_weights) =
      if c_local_error_set VS.! 0 == 0
        then (mempty, mempty)
        else coerce (c_local_error, c_var_weight)
    diagnostics = SundialsDiagnostics
      (fromIntegral $ c_diagnostics VS.!0)
      (fromIntegral $ c_diagnostics VS.!1)
      (fromIntegral $ c_diagnostics VS.!2)
      (fromIntegral $ c_diagnostics VS.!3)
      (fromIntegral $ c_diagnostics VS.!4)
      (fromIntegral $ c_diagnostics VS.!5)
      (fromIntegral $ c_diagnostics VS.!6)
      (fromIntegral $ c_diagnostics VS.!7)
      (fromIntegral $ c_diagnostics VS.!8)
      (fromIntegral $ c_diagnostics VS.!9)
      (toEnum . fromIntegral $ c_diagnostics VS.! 10)
  return $
    if ret == T.cV_SUCCESS
      then
        Right $ SundialsSolution
          { actualTimeGrid = extractTimeGrid output_mat
          , solutionMatrix = dropTimeGrid output_mat
          , diagnostics
          }
      else
        Left ErrorDiagnostics
          { partialResults = output_mat
          , errorCode = fromIntegral ret
          , errorEstimates = local_errors
          , varWeights = var_weights
          }
  where
    -- The time grid is the first column of the result matrix
    extractTimeGrid :: Matrix Double -> VS.Vector Double
    extractTimeGrid = head . toColumns
    dropTimeGrid :: Matrix Double -> Matrix Double
    dropTimeGrid = fromColumns . tail . toColumns

-- | The common solving logic between ARKode and CVode
solveCommon
  :: (IsMethod method, Katip m)
  => (CConsts -> CVars (VS.MVector RealWorld) -> LogEnv -> IO CInt)
      -- ^ the CVode/ARKode solving function; mostly inline-C code
  -> ODEOpts method
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveCommon solve_c opts problem@(OdeProblem{..})

  | VS.null odeInitCond = -- 0-dimensional (empty) system

    return . Right $ SundialsSolution
      { actualTimeGrid = odeSolTimes
      , solutionMatrix = (VS.length odeSolTimes >< 0) []
      , diagnostics = mempty
      }

  | otherwise = do

    log_env <- getLogEnv
    liftIO $ do -- the rest is in the IO monad
    vars <- allocateCVars problem
    ret <- withCConsts opts problem $ \consts ->
      solve_c consts vars log_env
    frozenVars <- freezeCVars vars
    assembleSolverResult problem ret frozenVars

-- | An auxiliary function to construct a storable vector from a C pointer
-- and length.
--
-- There doesn't seem to be a safe version of 'VS.unsafeFromForeignPtr0',
-- nor a way to clone an immutable vector, so we emulate it via an
-- intermediate mutable vector.
vecFromPtr
  :: VS.Storable a
  => Ptr a
  -> Int
  -> IO (VS.Vector a)
vecFromPtr ptr n = do
  fptr <- newForeignPtr_ ptr
  let mv = VSM.unsafeFromForeignPtr0 fptr n
  VS.freeze mv -- this does the copying and makes the whole thing safe

----------------------------------------------------------------------
--                           Logging
----------------------------------------------------------------------

-- | The Katip payload for logging Sundials errors
data SundialsErrorContext = SundialsErrorContext
  { sundialsErrorCode :: !Int
  , sundialsErrorModule :: !T.Text
  , sundialsErrorFunction :: !T.Text
  } deriving Generic
instance ToJSON SundialsErrorContext
instance ToObject SundialsErrorContext
instance LogItem SundialsErrorContext where
  payloadKeys _ _ = AllKeys

type ReportErrorFn =
  (  CInt    -- error code
  -> CString -- module name
  -> CString -- function name
  -> CString -- the message
  -> Ptr ()  -- user data (ignored)
  -> IO ()
  )

cstringToText :: CString -> IO T.Text
cstringToText = fmap T.decodeUtf8 . BS.packCString

reportErrorWithKatip :: LogEnv -> ReportErrorFn
reportErrorWithKatip log_env err_code c_mod_name c_func_name c_msg _userdata = do
  -- See Note [CV_TOO_CLOSE]
  if err_code == T.cV_TOO_CLOSE then pure () else do
  let
  mod_name <- cstringToText c_mod_name
  func_name <- cstringToText c_func_name
  msg <- cstringToText c_msg
  let
    severity :: Severity
    severity =
      if err_code <= 0
        then ErrorS
        else InfoS
    errCtx :: SundialsErrorContext
    errCtx = SundialsErrorContext
      { sundialsErrorCode = fromIntegral err_code
      , sundialsErrorModule = mod_name
      , sundialsErrorFunction = func_name
      }
  flip runReaderT log_env . unKatipT $ do
    logF errCtx "sundials" severity (logStr msg)

debugMsgWithKatip :: LogEnv -> CString -> IO ()
debugMsgWithKatip log_env cstr = do
  text <- cstringToText cstr
  flip runReaderT log_env . unKatipT $ do
    logF () "hmatrix-sundials" DebugS (logStr text)

-- From the former Types module

data EventHandlerResult = EventHandlerResult
  { eventStopSolver :: !Bool
    -- ^ should we stop the solver after handling this event?
  , eventRecord :: !Bool
    -- ^ should we record the state before and after the event in the ODE
    -- solution?
  , eventNewState :: !(VS.Vector Double)
    -- ^ the new state after the event has been applied
  }

type EventHandler
  =  Double -- ^ time
  -> VS.Vector Double -- ^ values of the variables
  -> VS.Vector Int
    -- ^ Vector of triggered event indices.
    -- If the vector is empty, this is a time-based event.
  -> IO EventHandlerResult

data OdeProblem = OdeProblem
  { odeEventConditions :: EventConditions
    -- ^ The event conditions
  , odeEventDirections :: V.Vector CrossingDirection
    -- ^ The requested directions of 0 crossing for each event. Also, the
    -- length of this vector tells us the number of events (even when
    -- 'odeEventConditions' is represented by a single C function).
  , odeMaxEvents :: !Int
    -- ^ The maximal number of events that may occur. This is needed to
    -- allocate enough space to store the events. If more events occur, an
    -- error is returned.
  , odeEventHandler :: EventHandler -- ^ The event handler.
  , odeTimeBasedEvents :: TimeEventSpec
  , odeRhs :: OdeRhs
    -- ^ The right-hand side of the system: either a Haskell function or
    -- a pointer to a compiled function.
  , odeJacobian :: Maybe OdeJacobian
    -- ^ The optional Jacobian (the arguments are the time and the state
    -- vector).
  , odeInitCond :: VS.Vector Double
    -- ^ The initial conditions of the problem.
  , odeSolTimes :: VS.Vector Double
    -- ^ The requested solution times. The actual solution times may be
    -- larger if any events occurred.
  , odeTolerances :: Tolerances
    -- ^ How much error is tolerated in each variable.
  }

data Tolerances = Tolerances
  { relTolerance :: !CDouble
  , absTolerances :: Either CDouble (VS.Vector CDouble)
    -- ^ If 'Left', then the same tolerance is used for all variables.
    --
    -- If 'Right', the vector should contain one tolerance per variable.
  } deriving (Show, Eq, Ord)

-- | The type of the C ODE RHS function.
type OdeRhsCType = CDouble -> Ptr SunVector -> Ptr SunVector -> Ptr UserData -> IO CInt

data UserData

-- | The right-hand side of an ODE system.
--
-- Can be either a Haskell function or a pointer to a C function.
data OdeRhs
  = OdeRhsHaskell (CDouble -> VS.Vector CDouble -> IO (VS.Vector CDouble))
  | OdeRhsC (FunPtr OdeRhsCType) (Ptr UserData)

-- | A version of 'OdeRhsHaskell' that accepts a pure function
odeRhsPure
  :: (CDouble -> VS.Vector CDouble -> VS.Vector CDouble)
  -> OdeRhs
odeRhsPure f = OdeRhsHaskell $ \t y -> return $ f t y

type OdeJacobianCType
  =  SunRealType   -- ^ @realtype t@
  -> Ptr SunVector -- ^ @N_Vector y@
  -> Ptr SunVector -- ^ @N_Vector fy@
  -> Ptr SunMatrix -- ^ @SUNMatrix Jac@
  -> Ptr UserData  -- ^ @void *user_data@
  -> Ptr SunVector -- ^ @N_Vector tmp1@
  -> Ptr SunVector -- ^ @N_Vector tmp2@
  -> Ptr SunVector -- ^ @N_Vector tmp3@
  -> IO CInt       -- ^ return value (0 if successful, >0 for a recoverable error, <0 for an unrecoverable error)

-- | The Jacobian of the right-hand side of an ODE system.
--
-- Can be either a Haskell function or a pointer to a C function.
data OdeJacobian
  = OdeJacobianHaskell (Double -> VS.Vector Double -> Matrix Double)
  | OdeJacobianC (FunPtr OdeJacobianCType)

data JacobianRepr
  = SparseJacobian !SparsePattern -- ^ sparse Jacobian with the given sparse pattern
  | DenseJacobian
  deriving (Show)

type EventConditionCType
  =  SunRealType     -- ^ @realtype t@
  -> Ptr SunVector   -- ^ @N_Vector y@
  -> Ptr SunRealType -- ^ @realtype *gout@
  -> Ptr UserData    -- ^ @void *user_data@
  -> IO CInt

data EventConditions
  = EventConditionsHaskell (Double -> VS.Vector Double -> VS.Vector Double)
  | EventConditionsC (FunPtr EventConditionCType)

-- | A way to construct 'EventConditionsHaskell' when there is no shared
-- computation among different functions
eventConditionsPure :: V.Vector (Double -> VS.Vector Double -> Double) -> EventConditions
eventConditionsPure conds = EventConditionsHaskell $ \t y ->
  V.convert $ V.map (\cond -> cond t y) conds

data ODEOpts method = ODEOpts {
    maxNumSteps :: Int32
  , minStep     :: Double
  , fixedStep   :: Double
      -- ^ If this is greater than 0.0, then a fixed-size step is used.
      --
      -- This is only recommended for testing/debugging, not for production
      -- use.
      --
      -- Also, this only has effect for ARKode; using this with CVode will
      -- trigger an error.
  , maxFail     :: Int32
  , odeMethod   :: method
  , initStep    :: Maybe Double
    -- ^ initial step size - by default, CVode
    -- estimates the initial step size to be the
    -- solution \(h\) of the equation
    -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
    -- \(\ddot{y}\) is an estimated value of the second
    -- derivative of the solution at \(t_0\)
  , jacobianRepr :: JacobianRepr
    -- ^ use a sparse matrix to represent the Jacobian
    -- and a sparse linear solver for Newton iterations
  } deriving (Show)

data SundialsDiagnostics = SundialsDiagnostics {
    odeGetNumSteps               :: Int
  , odeGetNumStepAttempts        :: Int
  , odeGetNumRhsEvals_fe         :: Int
  , odeGetNumRhsEvals_fi         :: Int
  , odeGetNumLinSolvSetups       :: Int
  , odeGetNumErrTestFails        :: Int
  , odeGetNumNonlinSolvIters     :: Int
  , odeGetNumNonlinSolvConvFails :: Int
  , dlsGetNumJacEvals            :: Int
  , dlsGetNumRhsEvals            :: Int
  , odeMaxEventsReached          :: Bool
  } deriving (Eq, Show, Generic, NFData)

instance Semigroup SundialsDiagnostics where
   (<>) (SundialsDiagnostics
          numSteps_1
          numStepAttempts_1
          numRhsEvals_fe_1
          numRhsEvals_fi_1
          numLinSolvSetups_1
          numErrTestFails_1
          numNonlinSolvIters_1
          numNonlinSolvConvFails_1
          numJacEvals_1
          numRhsEvals_1
          reachedMaxEvents_1)

        (SundialsDiagnostics
          numSteps_2
          numStepAttempts_2
          numRhsEvals_fe_2
          numRhsEvals_fi_2
          numLinSolvSetups_2
          numErrTestFails_2
          numNonlinSolvIters_2
          numNonlinSolvConvFails_2
          numJacEvals_2
          numRhsEvals_2
          reachedMaxEvents_2)

      = SundialsDiagnostics
          (numSteps_2 + numSteps_1)
          (numStepAttempts_2 + numStepAttempts_1)
          (numRhsEvals_fe_2 + numRhsEvals_fe_1)
          (numRhsEvals_fi_2 + numRhsEvals_fi_1)
          (numLinSolvSetups_2 + numLinSolvSetups_1)
          (numErrTestFails_2 + numErrTestFails_1)
          (numNonlinSolvIters_2 + numNonlinSolvIters_1)
          (numNonlinSolvConvFails_2 + numNonlinSolvConvFails_1)
          (numJacEvals_2 + numJacEvals_1)
          (numRhsEvals_2 + numRhsEvals_1)
          (reachedMaxEvents_1 || reachedMaxEvents_2)

instance Monoid SundialsDiagnostics
  where
    mempty = SundialsDiagnostics 0 0 0 0 0 0 0 0 0 0 False

data SundialsSolution =
  SundialsSolution
  { actualTimeGrid :: VS.Vector Double    -- ^ actual time grid returned by the solver (with duplicated event times)
  , solutionMatrix :: Matrix Double       -- ^ matrix of solutions: each column is an unknwown
  , diagnostics    :: SundialsDiagnostics -- ^ usual Sundials diagnostics
  }

data ErrorDiagnostics = ErrorDiagnostics
  { errorCode :: !Int
    -- ^ The numeric error code. Mostly useless at this point, since it is
    -- set to 1 under most error conditions. See 'solveOdeC'.
  , errorEstimates :: !(VS.Vector Double)
    -- ^ The local error estimates as returned by @CVodeGetEstLocalErrors@.
    -- Either an empty vector, or has the same dimensionality as the state
    -- space.
  , varWeights :: !(VS.Vector Double)
    -- ^ The weights with which errors are combined, equal to @1 / (atol_i + y_i * rtol)@.
    -- Either an empty vector, or has the same dimensionality as the state
    -- space.
  , partialResults :: !(Matrix Double)
    -- ^ Partial solution of the ODE system, up until the moment when
    -- solving failed. Contains the time as its first column.
  } deriving Show

-- | The direction in which a function should cross the x axis
data CrossingDirection = Upwards | Downwards | AnyDirection
  deriving (Generic, Eq, Show, NFData)

-- | A time-based event, implemented as an action that returns the time of
-- the next time-based event.
--
-- If there's an additional condition attached to a time-based event, it
-- should be verified in the event handler.
--
-- The action is supposed to be stateful, and the state of the action
-- should be updated by the event handler so that after a given time-based
-- event is handled, the action starts returning the time of the next
-- unhandled time-based event.
--
-- If there is no next time-based event, the action should return +Inf.
newtype TimeEventSpec = TimeEventSpec (IO Double)

sunTypesTable :: Map.Map TypeSpecifier TH.TypeQ
sunTypesTable = Map.fromList
  [
    (TypeName "sunindextype", [t| SunIndexType |] )
  , (TypeName "realtype",     [t| SunRealType |] )
  , (TypeName "N_Vector",     [t| Ptr SunVector |] )
  , (TypeName "SUNMatrix",    [t| Ptr SunMatrix |] )
  , (TypeName "UserData",     [t| UserData |] )
  ]

-- | Allows to map between Haskell and C types
sunCtx :: Context
sunCtx = mempty {ctxTypesTable = sunTypesTable}
