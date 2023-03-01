module Numeric.Sundials
  ( -- * The solving function
    solve
    -- * Problem specification
  , OdeRhsCType
  , OdeRhs(..)
  , odeRhsPure
  , OdeJacobianCType
  , OdeJacobian(..)
  , UserData
  , JacobianRepr(..)
  , SparsePattern(..)
  , OdeProblem(..)
    -- * Events
  , EventHandler
  , EventHandlerResult(..)
  , EventConditionCType
  , EventConditions(..)
  , eventConditionsPure
  , CrossingDirection(..)
  , TimeEventSpec(..)
    -- * Solution
  , SundialsDiagnostics(..)
  , ErrorDiagnostics(..)
  , SundialsSolution(..)
    -- * Solving options
  , ODEOpts(..)
  , OdeMethod(..)
  , ARK.ARKMethod(..)
  , CV.CVMethod(..)
  , allOdeMethods
  , IsMethod(..)
  , MethodType(..)
  , Tolerances(..)
    -- * Low-level types from sundials
  , SunVector(..)
  , SunMatrix(..)
  , SunIndexType
  , SunRealType
    -- * Offsets
    -- ** NVector
  , nvectorContentOffset
    -- ** NVector_SERIAL
  , nvectorContentSerialLengthOffset
  , nvectorContentSerialDataOffset
    -- ** SUNMatrix
  , sunmatrixContentOffset
    -- ** SUNMatrix_DENSE
  , sunmatrixContentDenseDataOffset
    -- ** SUNMatrix_SPARSE
  , sunmatrixContentSparseIndexvalsOffset
  , sunmatrixContentSparseIndexptrsOffset
  , sunmatrixContentSparseDataOffset
  , sunmatrixContentSparseNnzOffset
  ) where

import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import qualified Numeric.Sundials.CVode as CV
import qualified Numeric.Sundials.ARKode as ARK
import Katip

import Foreign.C.Types
import Foreign.Ptr
import Foreign.Storable (peek, poke)
import Numeric.Sundials.Foreign as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import Data.Maybe
import Data.Int
import Numeric.LinearAlgebra.HMatrix as H hiding (Vector)
import GHC.Prim
import GHC.Generics
import Control.Monad.IO.Class
import Control.Monad.Cont
import Control.Exception
import Control.DeepSeq

-- | A supported ODE solving method, either by CVode or ARKode
data OdeMethod
  = CVMethod CV.CVMethod
  | ARKMethod ARK.ARKMethod
  deriving (Eq, Ord, Show, Read, Generic)

-- | List of all supported ODE methods
allOdeMethods :: [OdeMethod]
allOdeMethods =
  (CVMethod <$> [minBound .. maxBound]) ++
  (ARKMethod <$> [minBound .. maxBound])

instance IsMethod OdeMethod where
  methodToInt = \case
    ARKMethod m -> methodToInt m
    CVMethod m -> methodToInt m
  methodType = \case
    ARKMethod m -> methodType m
    CVMethod m -> methodType m

data ODEOpts = ODEOpts {
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
  , odeMethod   :: OdeMethod
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

-- | Solve an ODE system using either ARKode or CVode (depending on what
-- @method@ is instantiated with).
solve
  :: forall m . Katip m
  => ODEOpts -- ^ solver options
  -> OdeProblem -- ^ the ODE system to solve
  -> m (Either ErrorDiagnostics SundialsSolution)
solve opts =
  let
    solveC =
      case odeMethod opts of
        CVMethod{} -> CV.solveC
        ARKMethod{} -> ARK.solveC
  in
    solveCommon solveC opts

-- | The common solving logic between ARKode and CVode
solveCommon
  :: (Katip m)
  => (CConsts -> CVars (VS.MVector RealWorld) -> LogEnv -> IO CInt)
      -- ^ the CVode/ARKode solving function; mostly inline-C code
  -> ODEOpts
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

withCConsts
  :: ODEOpts
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

              -- Force the evaluation and catch exception
              -- That's important because the following function will be called by sundials (C context).
              -- In this context, the RTS will crash in the event of an exception.
              -- NOTE: we force to NF using "force" in order to get all hidden
              -- exception. However, here, that's Storable Vector, which are
              -- strict on their argument, so forcing to NF is an O(1)
              -- operation, it only needs to force to WHNF.
              resM <- try (evaluate (Control.DeepSeq.force $ coerce f t y))
              case resM of
                Right res -> do
                  -- FIXME: We should be able to use poke somehow
                  T.vectorToC res (fromIntegral c_n_event_specs) out_ptr
                  return 0
                Left (_ :: SomeException) -> do
                  -- There was an exception, just return a non 0 value, so
                  -- sundial know that it must terminate the solving.
                  return 1
      funptr <- ContT $ bracket (mkEventConditionsC funIO) freeHaskellFunPtr
      return funptr
  return CConsts{..}
