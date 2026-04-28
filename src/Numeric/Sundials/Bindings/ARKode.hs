{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Numeric.Sundials.Bindings.ARKode where

import Data.Void
import Foreign
import Foreign.C
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import GHC.OverloadedLabels

-- | An opaque pointer to a ARKodeMem
newtype ARKodeMem = ARKodeMem (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "ARKStepCreate" cARKStepCreate :: FunPtr OdeRhsCType -> FunPtr OdeRhsCType -> CDouble -> N_Vector -> SUNContext -> IO ARKodeMem

foreign import ccall "ARKStepFree" cARKStepFree :: Ptr ARKodeMem -> IO ()

foreign import ccall "ARKodeSetUserData" cARKodeSetUserData :: ARKodeMem -> Ptr UserData -> IO (Flag ARKode)

foreign import ccall "ARKodeSetMinStep" cARKodeSetMinStep :: ARKodeMem -> CDouble -> IO (Flag ARKode)

foreign import ccall "ARKodeSetMaxStep" cARKodeSetMaxStep :: ARKodeMem -> CDouble -> IO (Flag ARKode)

foreign import ccall "ARKodeSetMaxNumSteps" cARKodeSetMaxNumSteps :: ARKodeMem -> SunIndexType -> IO (Flag ARKode)

foreign import ccall "ARKodeSetMaxErrTestFails" cARKodeSetMaxErrTestFails :: ARKodeMem -> CInt -> IO (Flag ARKode)

foreign import ccall "ARKodeSVtolerances" cARKodeSVtolerances :: ARKodeMem -> CDouble -> N_Vector -> IO (Flag ARKode)

foreign import ccall "ARKodeRootInit" cARKodeRootInit :: ARKodeMem -> CInt -> FunPtr EventConditionCType -> IO (Flag ARKode)

foreign import ccall "ARKodeSetRootDirection" cARKodeSetRootDirection :: ARKodeMem -> Ptr CInt -> IO (Flag ARKode)

foreign import ccall "ARKodeSetNoInactiveRootWarn" cARKodeSetNoInactiveRootWarn :: ARKodeMem -> IO (Flag ARKode)

foreign import ccall "ARKodeSetLinearSolver" cARKodeSetLinearSolver :: ARKodeMem -> SUNLinearSolver -> SUNMatrix -> IO (Flag ARKode)

foreign import ccall "ARKodeSetInitStep" cARKodeSetInitStep :: ARKodeMem -> CDouble -> IO (Flag ARKode)

foreign import ccall "ARKodeEvolve" cARKodeEvolve :: ARKodeMem -> CDouble -> N_Vector -> Ptr CDouble -> CInt -> IO (Flag ARKode)

foreign import ccall "ARKStepReInit"
  cARKStepReInit ::
    ARKodeMem ->
    FunPtr OdeRhsCType ->
    FunPtr a2 ->
    CDouble ->
    N_Vector ->
    IO (Flag ARKode)

foreign import ccall "ARKodeGetRootInfo" cARKodeGetRootInfo :: ARKodeMem -> Ptr CInt -> IO (Flag ARKode)

foreign import ccall "ARKodeSetJacFn" cARKodeSetJacFn :: ARKodeMem -> FunPtr OdeJacobianCType -> IO (Flag ARKode)

foreign import ccall "ARKodeSetFixedStep" cARKodeSetFixedStep :: ARKodeMem -> CDouble -> IO ()

-- Note: the CInt are actually Enum and this could be enforced
foreign import ccall "ARKStepSetTableNum" cARKStepSetTableNum :: ARKodeMem -> CInt -> CInt -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumSteps" cARKodeGetNumSteps :: ARKodeMem -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumLinSolvSetups" cARKodeGetNumLinSolvSetups :: ARKodeMem -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumErrTestFails" cARKodeGetNumErrTestFails :: ARKodeMem -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumNonlinSolvIters" cARKodeGetNumNonlinSolvIters :: ARKodeMem -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumNonlinSolvConvFails" cARKodeGetNumNonlinSolvConvFails :: ARKodeMem -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumJacEvals" cARKodeGetNumJacEvals :: ARKodeMem -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKStepGetNumRhsEvals" cARKStepGetNumRhsEvals :: ARKodeMem -> Ptr CLong -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumLinRhsEvals" cARKodeGetNumLinRhsEvals :: ARKodeMem -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumStepAttempts" cARKodeGetNumStepAttempts :: ARKodeMem -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumGEvals" cARKodeGetNumGEvals :: ARKodeMem -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetEstLocalErrors" cARKodeGetEstLocalErrors :: ARKodeMem -> N_Vector -> IO (Flag ARKode)

foreign import ccall "ARKodeGetErrWeights" cARKodeGetErrWeights :: ARKodeMem -> N_Vector -> IO (Flag ARKode)

data ARKode

instance IsLabel "SUCCESS" (Flag ARKode) where
  fromLabel = Flag ARK_SUCCESS

instance IsLabel "TOO_CLOSE" (Flag ARKode) where
  fromLabel = Flag ARK_TOO_CLOSE

instance IsLabel "ROOT_RETURN" (Flag ARKode) where
  fromLabel = Flag ARK_ROOT_RETURN
