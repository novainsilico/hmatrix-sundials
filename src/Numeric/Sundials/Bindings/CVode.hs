{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Numeric.Sundials.Bindings.CVode where

import Data.Void
import Foreign
import Foreign.C
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign

-- | An opaque pointer to a CVodeMem
newtype CVodeMem = CVodeMem (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "CVodeCreate" cCVodeCreate :: CInt -> SUNContext -> IO CVodeMem

foreign import ccall "CVodeFree" cCVodeFree :: Ptr CVodeMem -> IO ()

foreign import ccall "CVodeInit" cCVodeInit :: CVodeMem -> FunPtr OdeRhsCType -> SunRealType -> N_Vector -> IO CInt

foreign import ccall "CVodeSetUserData" cCVodeSetUserData :: CVodeMem -> Ptr UserData -> IO CInt

foreign import ccall "CVodeSetMinStep" cCVodeSetMinStep :: CVodeMem -> CDouble -> IO CInt

foreign import ccall "CVodeSetMaxStep" cCVodeSetMaxStep :: CVodeMem -> CDouble -> IO CInt

foreign import ccall "CVodeSetMaxNumSteps" cCVodeSetMaxNumSteps :: CVodeMem -> SunIndexType -> IO CInt

foreign import ccall "CVodeSetMaxErrTestFails" cCVodeSetMaxErrTestFails :: CVodeMem -> CInt -> IO CInt

foreign import ccall "CVodeSVtolerances" cCVodeSVtolerances :: CVodeMem -> CDouble -> N_Vector -> IO CInt

foreign import ccall "CVodeRootInit" cCVodeRootInit :: CVodeMem -> CInt -> FunPtr EventConditionCType -> IO CInt

foreign import ccall "CVodeSetRootDirection" cCVodeSetRootDirection :: CVodeMem -> Ptr CInt -> IO CInt

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

foreign import ccall "CVodeGetNumGEvals" cCVodeGetNumGEvals :: CVodeMem -> Ptr CLong -> IO CInt

foreign import ccall "CVodeGetEstLocalErrors" cCVodeGetEstLocalErrors :: CVodeMem -> N_Vector -> IO CInt

foreign import ccall "CVodeGetErrWeights" cCVodeGetErrWeights :: CVodeMem -> N_Vector -> IO CInt
