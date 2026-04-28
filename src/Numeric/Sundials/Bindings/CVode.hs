{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Numeric.Sundials.Bindings.CVode where

import Data.Void
import Foreign
import Foreign.C
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import GHC.OverloadedLabels

-- | An opaque pointer to a CVodeMem
newtype CVodeMem = CVodeMem (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "CVodeCreate" cCVodeCreate :: CInt -> SUNContext -> IO CVodeMem

foreign import ccall "CVodeFree" cCVodeFree :: Ptr CVodeMem -> IO ()

foreign import ccall "CVodeInit" cCVodeInit :: CVodeMem -> FunPtr OdeRhsCType -> SunRealType -> N_Vector -> IO (Flag CVode)

foreign import ccall "CVodeSetUserData" cCVodeSetUserData :: CVodeMem -> Ptr UserData -> IO (Flag CVode)

foreign import ccall "CVodeSetMinStep" cCVodeSetMinStep :: CVodeMem -> CDouble -> IO (Flag CVode)

foreign import ccall "CVodeSetMaxStep" cCVodeSetMaxStep :: CVodeMem -> CDouble -> IO (Flag CVode)

foreign import ccall "CVodeSetMaxNumSteps" cCVodeSetMaxNumSteps :: CVodeMem -> SunIndexType -> IO (Flag CVode)

foreign import ccall "CVodeSetMaxErrTestFails" cCVodeSetMaxErrTestFails :: CVodeMem -> CInt -> IO (Flag CVode)

foreign import ccall "CVodeSVtolerances" cCVodeSVtolerances :: CVodeMem -> CDouble -> N_Vector -> IO (Flag CVode)

foreign import ccall "CVodeRootInit" cCVodeRootInit :: CVodeMem -> CInt -> FunPtr EventConditionCType -> IO (Flag CVode)

foreign import ccall "CVodeSetRootDirection" cCVodeSetRootDirection :: CVodeMem -> Ptr CInt -> IO (Flag CVode)

foreign import ccall "CVodeSetNoInactiveRootWarn" cCVodeSetNoInactiveRootWarn :: CVodeMem -> IO (Flag CVode)

foreign import ccall "CVodeSetLinearSolver" cCVodeSetLinearSolver :: CVodeMem -> SUNLinearSolver -> SUNMatrix -> IO (Flag CVode)

foreign import ccall "CVodeSetInitStep" cCVodeSetInitStep :: CVodeMem -> CDouble -> IO (Flag CVode)

foreign import ccall "CVode" cCVode :: CVodeMem -> CDouble -> N_Vector -> Ptr CDouble -> CInt -> IO (Flag CVode)

foreign import ccall "CVodeReInit" cCVodeReInit :: CVodeMem -> CDouble -> N_Vector -> IO ()

foreign import ccall "CVodeGetRootInfo" cCVodeGetRootInfo :: CVodeMem -> Ptr CInt -> IO (Flag CVode)

foreign import ccall "CVodeSetJacFn" cCVodeSetJacFn :: CVodeMem -> FunPtr OdeJacobianCType -> IO (Flag CVode)

foreign import ccall "CVodeGetNumSteps" cCVodeGetNumSteps :: CVodeMem -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumLinSolvSetups" cCVodeGetNumLinSolvSetups :: CVodeMem -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumErrTestFails" cCVodeGetNumErrTestFails :: CVodeMem -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumNonlinSolvIters" cCVodeGetNumNonlinSolvIters :: CVodeMem -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumNonlinSolvConvFails" cCVodeGetNumNonlinSolvConvFails :: CVodeMem -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumJacEvals" cCVodeGetNumJacEvals :: CVodeMem -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumRhsEvals" cCVodeGetNumRhsEvals :: CVodeMem -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumLinRhsEvals" cCVodeGetNumLinRhsEvals :: CVodeMem -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumGEvals" cCVodeGetNumGEvals :: CVodeMem -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetEstLocalErrors" cCVodeGetEstLocalErrors :: CVodeMem -> N_Vector -> IO (Flag CVode)

foreign import ccall "CVodeGetErrWeights" cCVodeGetErrWeights :: CVodeMem -> N_Vector -> IO (Flag CVode)

data CVode

instance IsLabel "SUCCESS" (Flag CVode) where
  fromLabel = Flag CV_SUCCESS

instance IsLabel "TOO_CLOSE" (Flag CVode) where
  fromLabel = Flag CV_TOO_CLOSE

instance IsLabel "ROOT_RETURN" (Flag CVode) where
  fromLabel = Flag CV_ROOT_RETURN
