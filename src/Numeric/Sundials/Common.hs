-- | Common infrastructure for CVode/ARKode
{-# OPTIONS_GHC -Wno-name-shadowing #-}
module Numeric.Sundials.Common where

import Foreign.C.Types
import Foreign.Ptr
import Foreign.C.String
import Numeric.Sundials.Foreign as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import Numeric.LinearAlgebra.HMatrix as H hiding (Vector)
import Control.DeepSeq
import Katip
import Data.Aeson
import qualified Data.Text as T
import qualified Data.Text.Encoding as T
import qualified Data.ByteString as BS
import Control.Monad.Reader
import GHC.Generics (Generic)
import Foreign.ForeignPtr

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

allocateCVars :: OdeProblem -> IO (CVars (VS.MVector VSM.RealWorld))
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
freezeCVars :: CVars (VS.MVector VSM.RealWorld) -> IO (CVars VS.Vector)
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
  , c_is_differential :: VS.Vector CDouble
  -- ^ For IDA: tells if a value is driven by a differential (==1.0) or algebraic (==0.0)
  , c_init_differentials :: VS.Vector CDouble
  -- ^ For IDA: initial value for the differentials. That's a "best guess"
  , c_rhs :: FunPtr OdeRhsCType
  -- ^ For ode solving, ode RHS function
  , c_ida_res :: FunPtr IDAResFn
  -- ^ For IDA: residual function

  , c_rhs_userdata :: Ptr UserData
  , c_rtol :: CDouble
  , c_atol :: VS.Vector CDouble
  , c_n_event_specs :: CInt
  , c_event_fn :: FunPtr EventConditionCType
  -- ^ Root solving for ode problem
  , c_event_fn_ida :: FunPtr IDARootFn
  -- ^ Root solving for ida problem

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
  -- ^ Jacobian for ode problem
  , c_jac_ida :: FunPtr IDALsJacFn
  -- ^ Jacobian for IDA problem

  , c_sparse_jac :: CInt
      -- ^ If 0, use a dense matrix.
      -- If non-0, use a sparse matrix with that number of non-zero
      -- elements.
  , c_requested_event_direction :: VS.Vector CInt
  -- TODO: this function is not called from C, so we can be smart and actually
  -- use Maybe, instead of sentinel value
  , c_next_time_event :: IO CDouble
  , c_max_events :: CInt
  , c_minstep :: CDouble
  , c_maxstep :: Maybe CDouble
  , c_fixedstep :: CDouble
  , c_max_n_steps :: SunIndexType
  , c_max_err_test_fails :: CInt
  , c_init_step_size_set :: CInt
  , c_init_step_size :: CDouble
  , c_ontimepoint :: TimePointHandler
  }

data MethodType = Explicit | Implicit
  deriving (Show, Eq)

class IsMethod method where
  methodToInt :: method -> CInt
  methodType :: method -> MethodType

matrixToSunMatrix :: Matrix Double -> T.SunMatrix
matrixToSunMatrix m = T.SunMatrix { T.rows = nr, T.cols = nc, T.vals = vs }
  where
    nr = fromIntegral $ H.rows m
    nc = fromIntegral $ H.cols m
    vs = VS.unsafeCoerceVector . VS.concat $ toColumns m

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
  mkIDAResFn :: IDAResFn -> IO (FunPtr IDAResFn)

foreign import ccall "wrapper"
  mkOdeJacobianC :: OdeJacobianCType -> IO (FunPtr OdeJacobianCType)
foreign import ccall "wrapper"
  mkIDALsJacFn :: IDALsJacFn -> IO (FunPtr IDALsJacFn)

foreign import ccall "wrapper"
  mkEventConditionsC :: EventConditionCType -> IO (FunPtr EventConditionCType)

foreign import ccall "wrapper"
  mkIDARootFn :: IDARootFn -> IO (FunPtr IDARootFn)

assembleSolverResult
  :: OdeProblem
  -> CInt
  -> CVars VS.Vector
  -> IO (Either ErrorDiagnostics SundialsSolution)
assembleSolverResult OdeProblem{..} ret CVars{..} = do
  let
    dim = VS.length odeInitCond
    n_rows = fromIntegral . VS.head $ c_n_rows
    output_mat = reshape (dim + 1) . VS.unsafeCoerceVector . subVector 0 ((dim + 1) * n_rows) $ c_output_mat
    (local_errors, var_weights) =
      if c_local_error_set VS.! 0 == 0
        then (mempty, mempty)
        else (VS.unsafeCoerceVector c_local_error, VS.unsafeCoerceVector c_var_weight)
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

type ReportErrorFnNew =
                 CInt -- Line no
                -> Ptr CChar -- function name
                -> Ptr CChar -- file
                -> Ptr CChar -- message
                -> CInt -- err code
                -> Ptr () -- user data
                -> Ptr () -- sundial context
                -> IO ()

wrapErrorNewApi :: ReportErrorFn -> ReportErrorFnNew
wrapErrorNewApi f _lineNumber functionName moduleName errorMessage errorCode userData _sundialContext = f errorCode moduleName functionName errorMessage userData

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

debugMsgWithKatip :: LogEnv -> String -> IO ()
debugMsgWithKatip log_env text = do
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

-- | This callback will be called when a timepoint is saved
-- Maybe in the future we'll use this as a stream provider, but for now it is
-- only used for debuging purpose.
-- Note that this is NOT required anymore to build a stream abstraction,
-- because the main loop is in haskell, so it can just stream the timepoint
-- directly.
type TimePointHandler
  =  CInt -- ^ timepoint index
  -> IO ()

-- | Represents the inner function of the system. The solver can solve
-- "residual problem" (including algebraic equations), such as @f(t, v, v') = 0@ with the 'IDAMethod'
-- solver and simple "ode" problems, such as @dv/dt = f(t, v)@, using
-- 'CVMethod', 'ARKMethod' and 'IDAMethod'.
--
-- Note that you can provide an 'OdeFunctions' while using 'IDAMethod', the
-- solver will wrap the functions for you.
data ProblemFunctions = 
      -- | When solving a "ode" problem, such as @dv/dt = f(t, v)@
       -- The right-hand side of the system: either a Haskell function or
       -- a pointer to a compiled function. That's the @f@ function.
      OdeProblemFunctions OdeRhs
      |
      -- | When solving a "residual" problem, such as @f(t, v, v') = 0@
      ResidualProblemFunctions ResidualFunctions

-- | When solving a "residual" problem, such as @f(t, v, v') = 0@
data ResidualFunctions = ResidualFunctions {
        odeResidual :: OdeResidual
      -- ^ The residual function. Is available, when solving with IDA, it will be
      -- used instead of the 'odeRhs'. That's the "f" function.
      , odeDifferentials :: (VS.Vector Double)
      -- ^ Only when solving with IDA, tells if the unknown is driven by an ode
      -- (==1.0) or by an algebraic equation (==0.0). Filled as all derivatives if
      -- not provided.
      , odeInitialDifferentials :: (VS.Vector Double)
      -- ^ Only when solving with IDA, provides the initial derivatievs of the
      -- unknown. Will be computed by the solver if not provided.
    }

data OdeProblem = OdeProblem
  { odeEventConditions :: EventConditions
    -- ^ The event conditions. Used either for "ode" or "residual" problem.
    -- TODO: exposes the event condition function for residual, it allows event
    -- to test @yp@.
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
  , odeFunctions :: ProblemFunctions
  , odeJacobian :: Maybe OdeJacobian
    -- ^ The optional Jacobian (the arguments are the time and the state
    -- vector). Used either for "ode" or "residual" problems.
    -- TODO: expose the Jacobian for residual problem, it is parametrized by
    -- the `yp` parameter.
  , odeInitCond :: VS.Vector Double
    -- ^ The initial conditions of the problem.
  , odeSolTimes :: VS.Vector Double
    -- ^ The requested solution times. The actual solution times may be
    -- larger if any events occurred.
  , odeTolerances :: Tolerances
    -- ^ How much error is tolerated in each variable.
  , odeOnTimePoint :: Maybe TimePointHandler
  -- ^ This is called everytime the solver stores a timepoint
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

-- | The residual function for IDA solving
type IDAResFn =
   -- t
   CDouble ->
   -- y
   Ptr SunVector ->
   -- yp
   Ptr SunVector ->
   -- F(t, y, yp) (out parameter)
   Ptr SunVector ->
   -- user data
   Ptr UserData ->
   -- | Return code
   IO CInt

data UserData

-- | The right-hand side of an ODE system.
--
-- Can be either a Haskell function or a pointer to a C function.
data OdeRhs
  = OdeRhsHaskell (CDouble -> VS.Vector CDouble -> IO (VS.Vector CDouble))
  | OdeRhsC (FunPtr OdeRhsCType) (Ptr UserData)

data OdeResidual
  = OdeResidualHaskell (CDouble -> VS.Vector CDouble -> VS.Vector CDouble -> IO (VS.Vector CDouble))
  | OdeResidualC (FunPtr IDAResFn) (Ptr UserData)

-- | A version of 'OdeRhsHaskell' that accepts a pure function
odeRhsPure
  :: (CDouble -> VS.Vector CDouble -> VS.Vector CDouble)
  -> ProblemFunctions
odeRhsPure f = OdeProblemFunctions $ OdeRhsHaskell $ \t y -> return $ f t y

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

type IDALsJacFn
  =  SunRealType   -- ^ @realtype t@
  ->  SunRealType   -- ^ @realtype cj@
  -> Ptr SunVector -- ^ @N_Vector y@
  -> Ptr SunVector -- ^ @N_Vector yp@
  -> Ptr SunVector -- ^ @N_Vector r@
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

type IDARootFn
  =  SunRealType     -- ^ @realtype t@
  -> Ptr SunVector   -- ^ @N_Vector y@
  -> Ptr SunVector   -- ^ @N_Vector yp@
  -> Ptr SunRealType -- ^ @realtype *gout@
  -> Ptr UserData    -- ^ @void *user_data@
  -> IO CInt

data EventConditions
  = EventConditionsHaskell (Double -> VS.Vector Double -> IO (VS.Vector Double))
  | EventConditionsC (FunPtr EventConditionCType)

-- | A way to construct 'EventConditionsHaskell' when there is no shared
-- computation among different functions
eventConditionsPure :: V.Vector (Double -> VS.Vector Double -> Double) -> EventConditions
eventConditionsPure conds = EventConditionsHaskell $ \t y ->
  pure $ V.convert $ V.map (\cond -> cond t y) conds

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
  deriving (Show)

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

{-
  After the inline-c refactoring, there are many todos / things which can be improved:

- We have a few function (event handling, conditions, next time step, ...,
  which are wrapped to C and that's not necessary anymore.
- The early exit / ptrStop logic is not required anymore
- Debuging can be moved to full katip with finer grained control
- All the c_n_rows logic, output_ind, input_ind, event_ind logic can just be
  removed and we could work on a stream / list of output elements.
- The diagnostics could directly be generated in the relevant struct, instead
  of being pushed in an opaque vector.
- A lot of "int" can be turned into "Bool"
- A few vector copy which are implemneted with loop may benefit from `memcpy`
- Many many pieces of the solving loop can be shared between the different
- solver, which will help us introducing new solver (such as IDA).
-}


-- ** Sundial bindings

