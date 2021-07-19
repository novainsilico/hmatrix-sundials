module Numeric.Sundials
  ( -- * The solving function
    solve
    -- * Types
  , ARKMethod(..)
  , CVMethod(..)
  , IsMethod(..)
  , MethodType(..)
  , Solver(..)
  , OdeProblem(..)
  , EventHandler
  , EventHandlerResult(..)
  , Tolerances(..)
  , OdeRhsCType
  , OdeRhs(..)
  , odeRhsPure
  , OdeJacobianCType
  , OdeJacobian(..)
  , UserData
  , JacobianRepr(..)
  , SparsePattern(..)
  , ODEOpts(..)
  , SundialsDiagnostics(..)
  , ErrorDiagnostics(..)
  , SundialsSolution(..)
  , CrossingDirection(..)
  , EventConditionCType
  , EventConditions(..)
  , eventConditionsPure
  , TimeEventSpec(..)
  , SunVector(..)
  , SunMatrix(..)
  , SunIndexType
  , SunRealType
  , sunCtx
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
import Numeric.Sundials.CVode as CV
import Numeric.Sundials.ARKode as ARK
import Katip

-- | Solve an ODE system using either ARKode or CVode (depending on what
-- @method@ is instantiated with).
solve
  :: forall m method . (Katip m, IsMethod method)
  => ODEOpts method -- ^ solver options
  -> OdeProblem -- ^ the ODE system to solve
  -> m (Either ErrorDiagnostics SundialsSolution)
solve =
  case methodSolver @method of
    CVode -> solveCommon CV.solveC
    ARKode -> solveCommon ARK.solveC
