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
  , IDA.IDAMethod(..)
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
import qualified Numeric.Sundials.IDA as IDA
import Katip

import Foreign.C.Types
import Foreign.Ptr
import Foreign.Storable (peek, poke)
import Numeric.Sundials.Foreign as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import Data.Maybe
import Data.Int
import Numeric.LinearAlgebra.HMatrix as H hiding ((<>), Vector)
import GHC.Generics
import Control.Monad.IO.Class
import Control.Monad.Cont
import Control.Exception
import Control.Concurrent.MVar
import Data.Vector.Mutable (RealWorld)
import Data.Coerce (coerce)
import Numeric.Sundials.Bindings.Sundials (withSUNContext, withNVector_Serial, cNV_Ith_S, N_Vector (..))

-- | A supported ODE solving method, either by CVode or ARKode
data OdeMethod
  = CVMethod CV.CVMethod
  | ARKMethod ARK.ARKMethod
  | IDAMethod IDA.IDAMethod
  deriving (Eq, Ord, Show, Read, Generic)

-- | List of all supported ODE methods
allOdeMethods :: [OdeMethod]
allOdeMethods =
  (CVMethod <$> [minBound .. maxBound]) ++
  (ARKMethod <$> [minBound .. maxBound]) ++ 
  (IDAMethod <$> [minBound .. maxBound])

instance IsMethod OdeMethod where
  methodToInt = \case
    ARKMethod m -> methodToInt m
    CVMethod m -> methodToInt m
    IDAMethod m -> methodToInt m
  methodType = \case
    ARKMethod m -> methodType m
    CVMethod m -> methodType m
    IDAMethod m -> methodType m

data ODEOpts = ODEOpts {
    maxNumSteps :: Int32
  , minStep     :: Double
  , maxStep     :: Maybe Double
  -- ^ an optional max step size, the default value is +∞
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
        IDAMethod{} -> IDA.solveC

  in solveCommon solveC opts

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

-- | Classify the method required as if it is ode or residual
data ProblemType =
  -- | An ode method only contains @dv/dt = f(t, v)@ odes
  Ode
  |
  -- | A residual method, such as IDA, can solve @f(t, v, v') = 0@
  Residual

-- | Classify the method required as if it is ode or residual
getProblemType :: OdeMethod -> ProblemType
getProblemType (IDAMethod _) = Residual
getProblemType (ARKMethod _) = Ode
getProblemType (CVMethod _) = Ode


withCConsts
  :: ODEOpts
  -> OdeProblem
  -> (CConsts -> IO r)
  -> IO r
withCConsts ODEOpts{..} OdeProblem{..} = runContT $ do
  exceptionRef <- ContT $ bracket newEmptyMVar $ \mvar -> do
    -- If an exception happened during the execution of the solver in a
    -- Haskell FFI call.
    -- a) It is caught by the haskell execution context
    -- b) It is put in the mvar
    -- c) The haskell execution context returns 1 to the sundial solver
    --    This is interpreted as an unrecoverable error
    -- d) Sundials terminate
    -- e) The exception is reread here and rethrown
    --
    -- Note that we do that for any kind of exception, sync or async, because
    -- in anycase, it is rethrown.
    res <- tryReadMVar mvar
    case res of
      Nothing -> pure ()
      Just e -> throwIO e

  let
    dim = VS.length (c_init_cond :: VS.Vector SunRealType)
    c_init_cond = VS.unsafeCoerceVector odeInitCond
    c_dim = fromIntegral dim
    c_n_sol_times = fromIntegral . VS.length $ odeSolTimes
    c_sol_time = VS.unsafeCoerceVector odeSolTimes
    c_rtol = relTolerance odeTolerances
    c_atol = either (VS.replicate dim) id $ absTolerances odeTolerances
    c_minstep = coerce minStep
    c_maxstep = coerce <$> maxStep
    c_fixedstep = coerce fixedStep
    c_max_n_steps = fromIntegral maxNumSteps
    c_max_err_test_fails = fromIntegral maxFail
    c_init_step_size_set = fromIntegral . fromEnum $ isJust initStep
    c_init_step_size =  maybe 0 coerce initStep
    c_n_event_specs = fromIntegral $ V.length odeEventDirections
    c_requested_event_direction = V.convert $ V.map directionToInt odeEventDirections
    -- TODO: this is not called from sundial, so the "C" wrapping is not mandatory
    c_apply_event n_events event_indices_ptr t y_ptr y'_ptr stop_solver_ptr record_event_ptr = do
      event_indices <- vecFromPtr event_indices_ptr (fromIntegral n_events)
      y_vec <- peek y_ptr
      
      saveExceptionContext exceptionRef $ do
        EventHandlerResult{..} <-
          odeEventHandler
            (coerce t :: Double)
            (VS.unsafeCoerceVector $ sunVecVals y_vec :: VS.Vector Double)
            (VS.map fromIntegral event_indices :: VS.Vector Int)
        poke y'_ptr $ SunVector
          { sunVecN = sunVecN y_vec
          , sunVecVals = VS.unsafeCoerceVector eventNewState
          }
        poke stop_solver_ptr . fromIntegral $ fromEnum eventStopSolver
        poke record_event_ptr . fromIntegral $ fromEnum eventRecord
    c_max_events = fromIntegral odeMaxEvents
    c_next_time_event = do
      resM <- saveExceptionContextM exceptionRef (evaluate =<< coerce odeTimeBasedEvents)
      case resM of
        Just t -> pure t
        -- -1 for next event will be interpereted by the solver as a stop
        Nothing -> pure (-1)

    c_jac_set = fromIntegral . fromEnum $ isJust odeJacobian
    c_sparse_jac = case jacobianRepr of
      SparseJacobian (T.SparsePattern spat) ->
        VS.sum (VS.map fromIntegral spat) +
        -- additionally, add diagonal zeros, as they'll be allocated too
        sum [ if spat VS.! (i + i * dim) == 0 then 1 else 0 | i <- [0 .. dim-1] ]
      DenseJacobian -> 0
    c_method = methodToInt odeMethod

  -- Starting here, the code generates the different 'FunPtr' required by sundials, namely:
  --
  --    - the ode or residual function: c_rhs, c_ida_res
  --    - the event function: c_event_fn, c_event_fn_ida
  --    - the jacobian function: c_jac, c_jac_id,
  --    - IDA only data: c_is_differential, c_init_differentials
  --
  -- The code is super convoluted because of the different uses cases we aim to support:
  --
  -- - For ode solving: c_rhs, c_event_fn, c_jac
  -- - For ida solving: c_ida_res, c_event_fn_ida, c_jac_ida, c_is_differential, c_init_differentials
  -- - for ida solving, but when on ODE problem is required, we generate:
  --    - c_ida_res, by wrapping c_rhs
  --    - c_is_differential is set to 1.0 for each components
  --    - c_init_differentials is initialised using c_rhs
  --
  -- The latest point is done in order to test the solver on "non-residual"
  -- problem.
  --
  -- Most of the not required values are set to dummy values. This is imperfect
  -- and unsafe, we just assume that the test coverage does correctly test
  -- theses cases.
  (c_rhs, c_ida_res, c_rhs_userdata, c_is_differential, c_init_differentials) <- case odeFunctions of
    OdeProblemFunctions odeRhs -> do
      -- Estimation of the initial differentials, because not provided by the user
      let compute_initial_differentials c_rhs c_rhs_userdata = do
           -- If we don't know, let's just evaluate the rhs once
           let t0 = fromMaybe (error "no t0") $ c_sol_time VS.!? 0
           withSUNContext $ \sunctx -> do
                withNVector_Serial c_dim sunctx 6471 $ \y0_ptr -> do
                withNVector_Serial c_dim sunctx 6471 $ \res_ptr -> do
                  VS.imapM_ (\i v -> cNV_Ith_S y0_ptr i v) c_init_cond
                  resC <- (runOdeRhs c_rhs) t0 (getSunVector y0_ptr) (getSunVector res_ptr) c_rhs_userdata
                  if resC /= 0
                  then
                    error $ "sundials returned an error during initialisation phase: " <> show resC
                  else do
                    res <- sunVecVals <$> peek (getSunVector res_ptr)
                    pure res

      case odeRhs of
        -- TODO: maybe we can leverage the user data somewhere in order to work
        -- with DDE and other stuffs and "constant" values (e.g. time-varying
        -- categoricals for example)
        OdeRhsC ptr u -> do
          case getProblemType odeMethod of
            Residual -> do
              -- If we don't know, let's assume that everything is differential
              let c_is_differential = VS.replicate dim 1.0
              funptrida <- wrap_ide_ode_rhs ptr
              initDifferentials <- liftIO $ compute_initial_differentials ptr u
              return (ptr, funptrida, u, c_is_differential, initDifferentials)
            Ode -> do
              return (ptr, nullFunPtr, u, mempty, mempty)
        OdeRhsHaskell fun -> do
          let
            funIO :: OdeRhsCType
            funIO t y f _ptr = do
              sv <- peek y

              -- Save the exception (if any)
              saveExceptionContext exceptionRef $ do
                r <- fun t (sunVecVals sv)

                -- Note: the following operation will force "r"
                -- and discover any hidden exception
                poke f $ SunVector { sunVecN = sunVecN sv
                                   , sunVecVals = r
                                   }

            -- In case the user does not provide a residual function, we build
            -- one from the ode rhs provided function.
            funIdaCompatIO :: IDAResFn
            funIdaCompatIO t y yp f _ptr = do
              -- Save the exception (if any)
              saveExceptionContext exceptionRef $ do
                sv <- peek y
                svp <- peek yp

                ypComputed <- fun t (sunVecVals sv)
                -- The residual function is F(y, yp, t) = 0
                -- However, we only have yp_rhs = f(y, t)
                --
                -- So we build F(y, yp, t) = yp_rhs - yp = f(y, t) - yp
                let res = VS.zipWith (-) ypComputed (sunVecVals svp)

                -- Note: the following operation will force "res"
                -- and discover any hidden exception
                poke f $ SunVector { sunVecN = sunVecN sv
                                   , sunVecVals = res
                                   }

          case getProblemType odeMethod of
            Residual -> do
              -- We solve an IDA problem, we don't care about the ode
              -- implementation, however the rhs can be useful when evaluating
              -- initial derivatives
              funptr <- ContT $ bracket (mkOdeRhsC funIO) freeHaskellFunPtr
              funidaptr <- ContT $ bracket (mkIDAResFn funIdaCompatIO) freeHaskellFunPtr

              -- If we don't know, let's assume that everything is differential
              let c_is_differential = VS.replicate dim 1.0
              initDifferentials <- liftIO $ compute_initial_differentials funptr nullPtr
              return (funptr, funidaptr, nullPtr, c_is_differential, initDifferentials)
            Ode -> do
              -- We will solve an ode problem, we don't care about the ida implementation
              funptr <- ContT $ bracket (mkOdeRhsC funIO) freeHaskellFunPtr
              let funidaptr = nullFunPtr
              -- We don't care about differential informations
              let c_is_differential = mempty
              return (funptr, funidaptr, nullPtr, c_is_differential, mempty)
    ResidualProblemFunctions ResidualFunctions{..} -> do
          (funidaptr, userdataptr) <- case odeResidual of
             OdeResidualHaskell odeResidualF -> do
              let
                -- That's a correct residual function
                funIdaResidualIO = fn
                  where 
                    fn t y yp residual _ptr = do
                       -- Save the exception (if any)
                       saveExceptionContext exceptionRef $ do
                         sv <- peek y
                         svp <- peek yp
                         res <- odeResidualF t (sunVecVals sv) (sunVecVals svp)

                         -- Note: the following operation will force "res"
                         -- and discover any hidden exception
                         poke residual $ SunVector { sunVecN = sunVecN sv
                                            , sunVecVals = res
                                            }
              funptr <- ContT $ bracket (mkIDAResFn funIdaResidualIO) freeHaskellFunPtr
              pure (funptr, nullPtr)
             OdeResidualC funptr userdataptr -> pure (funptr, userdataptr)
          return (nullFunPtr, funidaptr, userdataptr, VS.unsafeCoerceVector odeDifferentials, VS.unsafeCoerceVector odeInitialDifferentials)
  let c_ontimepoint idx = do
        case odeOnTimePoint of
          Nothing -> pure ()
          Just fun -> do
            -- Save the exception (if any)
            _ <- saveExceptionContext exceptionRef $ do
              fun idx
            pure ()
  (c_jac, c_jac_ida) <-
    case odeJacobian of
      Nothing   -> return (nullFunPtr, nullFunPtr)
      Just (OdeJacobianC fptr) -> do
        funptrida <- wrap_ide_ode_jacobian fptr
        return (fptr, funptrida)
      Just (OdeJacobianHaskell jac_fn) -> do
      let
        funIO :: OdeJacobianCType
        funIO t y_ptr _fy_ptr jac_ptr _userdata _tmp1 _tmp2 _tmp3 = do
          y <- peek y_ptr
          let jac = matrixToSunMatrix $
                jac_fn
                  (coerce t :: Double)
                  (VS.unsafeCoerceVector $ sunVecVals y :: VS.Vector Double)
          case jacobianRepr of
            DenseJacobian -> poke jac_ptr jac
            SparseJacobian spat -> poke (castPtr jac_ptr) (T.SparseMatrix spat jac)
          return 0
        funIdaIO :: IDALsJacFn
        funIdaIO t cj y_ptr _yp_ptr _r jac_ptr _userdata _tmp1 _tmp2 _tmp3 = do
          y <- peek y_ptr
          let jac = matrixToSunMatrix $
                (jac_fn
                  (coerce t :: Double)
                  (VS.unsafeCoerceVector $ sunVecVals y :: VS.Vector Double))
                     - (realToFrac cj * ident (VS.length (sunVecVals y)))
          case jacobianRepr of
            DenseJacobian -> poke jac_ptr jac
            SparseJacobian spat -> poke (castPtr jac_ptr) (T.SparseMatrix spat jac)
          return 0
      case getProblemType odeMethod of
        Ode -> do
          funptr <- ContT $ bracket (mkOdeJacobianC funIO) freeHaskellFunPtr
          let funidaptr = nullFunPtr
          return (funptr, funidaptr)
        Residual -> do
          let funptr = nullFunPtr
          funidaptr <- ContT $ bracket (mkIDALsJacFn funIdaIO) freeHaskellFunPtr
          return (funptr, funidaptr)

  (c_event_fn, c_event_fn_ida) <-
    case odeEventConditions of
      EventConditionsC fptr -> do
        case getProblemType odeMethod of
          Ode -> do
            let funptrida = nullFunPtr
            return (fptr, funptrida)
          Residual -> do
            funptrida <- wrap_ida_event_condition fptr
            return (nullFunPtr, funptrida)
      EventConditionsHaskell f -> do
        let
          funIO :: EventConditionCType
          funIO t y_ptr out_ptr _ptr = do
                y <- sunVecVals <$> peek y_ptr

                saveExceptionContext exceptionRef $ do
                   res <- f (coerce t) (VS.unsafeCoerceVector y)
                   -- FIXME: We should be able to use poke somehow
                   -- Note: the following operation will force "res"
                   -- and discover any hidden exception
                   T.vectorToC (VS.unsafeCoerceVector res) (fromIntegral c_n_event_specs) out_ptr

          -- TODO: the yp_ptr could be used in root functions
          funIdaIO :: IDARootFn
          funIdaIO t y_ptr _yp_ptr out_ptr _ptr = do
                y <- sunVecVals <$> peek y_ptr

                saveExceptionContext exceptionRef $ do
                   res <- f (coerce t) (VS.unsafeCoerceVector y)
                   -- FIXME: We should be able to use poke somehow
                   -- Note: the following operation will force "res"
                   -- and discover any hidden exception
                   T.vectorToC (VS.unsafeCoerceVector res) (fromIntegral c_n_event_specs) out_ptr
        case getProblemType odeMethod of
          Ode -> do
            funptr <- ContT $ bracket (mkEventConditionsC funIO) freeHaskellFunPtr
            let funidaptr = nullFunPtr
            return (funptr, funidaptr)
          Residual -> do
            let funptr = nullFunPtr
            funidaptr <- ContT $ bracket (mkIDARootFn funIdaIO) freeHaskellFunPtr
            return (funptr, funidaptr)
      EventConditionsResidualC fptr -> do
        case getProblemType odeMethod of
          Ode -> do
            error "Cannot solve a system without IDA when the event condition requires yp"
          Residual -> do
            return (nullFunPtr, fptr)
      EventConditionsResidualHaskell f -> do
      let
        funIdaIO :: IDARootFn
        funIdaIO t y_ptr yp_ptr out_ptr _ptr = do
              y <- sunVecVals <$> peek y_ptr
              yp <- sunVecVals <$> peek yp_ptr

              saveExceptionContext exceptionRef $ do
                 res <- f (coerce t) (VS.unsafeCoerceVector y) (VS.unsafeCoerceVector yp)
                 -- FIXME: We should be able to use poke somehow
                 -- Note: the following operation will force "res"
                 -- and discover any hidden exception
                 T.vectorToC (VS.unsafeCoerceVector res) (fromIntegral c_n_event_specs) out_ptr
      case getProblemType odeMethod of
        Ode -> do
          error "Cannot solve a system without IDA when the event condition requires yp"
        Residual -> do
          let funptr = nullFunPtr
          funidaptr <- ContT $ bracket (mkIDARootFn funIdaIO) freeHaskellFunPtr
          return (funptr, funidaptr)

  return CConsts{..}


-- | Wrapped to call the event condition directly from haskell code. This is
-- used to wrap the event condition "ode" style in a "residual" style.
foreign import ccall "dynamic"
  runOdeEventCondition:: FunPtr EventConditionCType  -> EventConditionCType 

-- | Convert event condition from "ode" to "residual" api.
wrap_ida_event_condition :: FunPtr EventConditionCType -> ContT r IO (FunPtr IDARootFn)
wrap_ida_event_condition funptr = do
  let 
        funIdaIO t y _yp res userdata = do
          (runOdeEventCondition funptr) t y res userdata
  funidaptr <- ContT $ bracket (mkIDARootFn funIdaIO) freeHaskellFunPtr
  pure funidaptr

-- | Wrapper to call the ode rhs directly from haskell code. This is used to
-- compute the initial differential in IDA solving when they are not provided
-- by the user.
foreign import ccall "dynamic"
  runOdeRhs :: FunPtr OdeRhsCType -> OdeRhsCType

-- | Convert ode rhs to residual function
wrap_ide_ode_rhs :: FunPtr OdeRhsCType -> ContT r IO (FunPtr IDAResFn)
wrap_ide_ode_rhs funptr = do
  let 
        funIdaIO t y yp res userdata = do
          -- The residual function is F(y, yp, t) = 0
          -- However, we only have yp_rhs = f(y, t)
          --
          -- So we build F(y, yp, t) = yp_rhs - yp = f(y, t) - yp
          ypVec <- peek yp
          error_code <- (runOdeRhs funptr) t y res userdata
          resVec <- peek res
          let res' = VS.zipWith (-) (sunVecVals resVec) (sunVecVals ypVec)

          poke res $ SunVector
             { sunVecN = sunVecN resVec
             , sunVecVals = VS.unsafeCoerceVector res'
           }
          pure error_code

  funidaptr <- ContT $ bracket (mkIDAResFn funIdaIO) freeHaskellFunPtr
  pure funidaptr

foreign import ccall "dynamic"
  runOdeJac :: FunPtr OdeJacobianCType -> OdeJacobianCType

-- | Convert ode style jacobian function to residual
wrap_ide_ode_jacobian :: FunPtr OdeJacobianCType -> ContT r IO (FunPtr IDALsJacFn)
wrap_ide_ode_jacobian funptr = do
  let 
        funIdaIO t _cj y_ptr _yp_ptr r jac_ptr _userdata _tmp1 _tmp2 _tmp3 = do
          (runOdeJac funptr) t y_ptr r jac_ptr _userdata _tmp1 _tmp2 _tmp3
  funidaptr <- ContT $ bracket (mkIDALsJacFn funIdaIO) freeHaskellFunPtr
  pure funidaptr


-- | Ensure that the @io@ does not raise any exception.
-- If it raises, save the exception in @exceptionRef@ and return an error code
-- to sundial, which will then terminate gracefully. The top level code is
-- responsible to raise the exception at the end of the solving.
saveExceptionContext :: MVar SomeException -> IO () -> IO CInt
saveExceptionContext exceptionRef io = do
  resM <- saveExceptionContextM exceptionRef io
  case resM of
    -- There was an exception, just return a non 0 value, so
    -- sundial know that it must terminate the solving.
    Nothing -> pure 1
    -- 0 is the "everything ok" code for sundial
    Just _ -> pure 0

saveExceptionContextM :: MVar SomeException -> IO a -> IO (Maybe a)
saveExceptionContextM exceptionRef io = do
  resM <- try io
  case resM of
    Right res -> do
      pure (Just res)
    Left (e :: SomeException) -> do
      -- Save the exception
      -- Note: it is possible that an exception was already saved, we just discard it then
      -- This seems to happen with IDA solver, when it fails, it will retry multiple times
      -- TODO: maybe we can compare the exception and store a deduplicated list with more context?
      _ <- tryPutMVar exceptionRef e
      pure Nothing
