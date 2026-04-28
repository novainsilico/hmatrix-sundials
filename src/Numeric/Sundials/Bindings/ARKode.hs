{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Numeric.Sundials.Bindings.ARKode where

import Data.Void
import Foreign
import Foreign.C
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign (SunIndexType)

-- | An opaque pointer to a ARKodeMem
newtype ARKodeMem = ARKodeMem (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "ARKStepCreate" cARKStepCreate :: FunPtr OdeRhsCType -> FunPtr OdeRhsCType -> CDouble -> N_Vector -> SUNContext -> IO ARKodeMem

foreign import ccall "ARKStepFree" cARKStepFree :: Ptr ARKodeMem -> IO ()

foreign import ccall "ARKodeSetUserData" cARKodeSetUserData :: ARKodeMem -> Ptr UserData -> IO CInt

foreign import ccall "ARKodeSetMinStep" cARKodeSetMinStep :: ARKodeMem -> CDouble -> IO CInt

foreign import ccall "ARKodeSetMaxStep" cARKodeSetMaxStep :: ARKodeMem -> CDouble -> IO CInt

foreign import ccall "ARKodeSetMaxNumSteps" cARKodeSetMaxNumSteps :: ARKodeMem -> SunIndexType -> IO CInt

foreign import ccall "ARKodeSetMaxErrTestFails" cARKodeSetMaxErrTestFails :: ARKodeMem -> CInt -> IO CInt

foreign import ccall "ARKodeSVtolerances" cARKodeSVtolerances :: ARKodeMem -> CDouble -> N_Vector -> IO CInt

foreign import ccall "ARKodeRootInit" cARKodeRootInit :: ARKodeMem -> CInt -> FunPtr EventConditionCType -> IO CInt

foreign import ccall "ARKodeSetRootDirection" cARKodeSetRootDirection :: ARKodeMem -> Ptr CInt -> IO CInt

foreign import ccall "ARKodeSetNoInactiveRootWarn" cARKodeSetNoInactiveRootWarn :: ARKodeMem -> IO CInt

foreign import ccall "ARKodeSetLinearSolver" cARKodeSetLinearSolver :: ARKodeMem -> SUNLinearSolver -> SUNMatrix -> IO CInt

foreign import ccall "ARKodeSetInitStep" cARKodeSetInitStep :: ARKodeMem -> CDouble -> IO CInt

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

foreign import ccall "ARKodeGetNumSteps" cARKodeGetNumSteps :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumLinSolvSetups" cARKodeGetNumLinSolvSetups :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumErrTestFails" cARKodeGetNumErrTestFails :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumNonlinSolvIters" cARKodeGetNumNonlinSolvIters :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumNonlinSolvConvFails" cARKodeGetNumNonlinSolvConvFails :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumJacEvals" cARKodeGetNumJacEvals :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKStepGetNumRhsEvals" cARKStepGetNumRhsEvals :: ARKodeMem -> Ptr CLong -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumLinRhsEvals" cARKodeGetNumLinRhsEvals :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumStepAttempts" cARKodeGetNumStepAttempts :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetNumGEvals" cARKodeGetNumGEvals :: ARKodeMem -> Ptr CLong -> IO CInt

foreign import ccall "ARKodeGetEstLocalErrors" cARKodeGetEstLocalErrors :: ARKodeMem -> N_Vector -> IO CInt

foreign import ccall "ARKodeGetErrWeights" cARKodeGetErrWeights :: ARKodeMem -> N_Vector -> IO CInt
