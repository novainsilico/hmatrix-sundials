{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Numeric.Sundials.Bindings.ARKode where

import Foreign
import Foreign.C
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import GHC.OverloadedLabels

foreign import ccall "ARKStepCreate" cARKStepCreate :: FunPtr OdeRhsCType -> FunPtr OdeRhsCType -> CDouble -> N_Vector -> SUNContext -> IO (SolverObject ARKode)

foreign import ccall "ARKStepFree" cARKStepFree :: Ptr (SolverObject ARKode) -> IO ()

foreign import ccall "ARKodeSetUserData" cARKodeSetUserData :: (SolverObject ARKode) -> Ptr UserData -> IO (Flag ARKode)

foreign import ccall "ARKodeSetMinStep" cARKodeSetMinStep :: (SolverObject ARKode) -> CDouble -> IO (Flag ARKode)

foreign import ccall "ARKodeSetMaxStep" cARKodeSetMaxStep :: (SolverObject ARKode) -> CDouble -> IO (Flag ARKode)

foreign import ccall "ARKodeSetMaxNumSteps" cARKodeSetMaxNumSteps :: (SolverObject ARKode) -> SunIndexType -> IO (Flag ARKode)

foreign import ccall "ARKodeSetMaxErrTestFails" cARKodeSetMaxErrTestFails :: (SolverObject ARKode) -> CInt -> IO (Flag ARKode)

foreign import ccall "ARKodeSVtolerances" cARKodeSVtolerances :: (SolverObject ARKode) -> CDouble -> N_Vector -> IO (Flag ARKode)

foreign import ccall "ARKodeRootInit" cARKodeRootInit :: (SolverObject ARKode) -> CInt -> FunPtr EventConditionCType -> IO (Flag ARKode)

foreign import ccall "ARKodeSetRootDirection" cARKodeSetRootDirection :: (SolverObject ARKode) -> Ptr CInt -> IO (Flag ARKode)

foreign import ccall "ARKodeSetNoInactiveRootWarn" cARKodeSetNoInactiveRootWarn :: (SolverObject ARKode) -> IO (Flag ARKode)

foreign import ccall "ARKodeSetLinearSolver" cARKodeSetLinearSolver :: (SolverObject ARKode) -> SUNLinearSolver -> SUNMatrix -> IO (Flag ARKode)

foreign import ccall "ARKodeSetInitStep" cARKodeSetInitStep :: (SolverObject ARKode) -> CDouble -> IO (Flag ARKode)

foreign import ccall "ARKodeEvolve" cARKodeEvolve :: (SolverObject ARKode) -> CDouble -> N_Vector -> Ptr CDouble -> CInt -> IO (Flag ARKode)

foreign import ccall "ARKStepReInit"
  cARKStepReInit ::
    (SolverObject ARKode) ->
    FunPtr OdeRhsCType ->
    FunPtr a2 ->
    CDouble ->
    N_Vector ->
    IO (Flag ARKode)

foreign import ccall "ARKodeGetRootInfo" cARKodeGetRootInfo :: (SolverObject ARKode) -> Ptr CInt -> IO (Flag ARKode)

foreign import ccall "ARKodeSetJacFn" cARKodeSetJacFn :: (SolverObject ARKode) -> FunPtr OdeJacobianCType -> IO (Flag ARKode)

foreign import ccall "ARKodeSetFixedStep" cARKodeSetFixedStep :: (SolverObject ARKode) -> CDouble -> IO ()

-- Note: the CInt are actually Enum and this could be enforced
foreign import ccall "ARKStepSetTableNum" cARKStepSetTableNum :: (SolverObject ARKode) -> CInt -> CInt -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumSteps" cARKodeGetNumSteps :: (SolverObject ARKode) -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumLinSolvSetups" cARKodeGetNumLinSolvSetups :: (SolverObject ARKode) -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumErrTestFails" cARKodeGetNumErrTestFails :: (SolverObject ARKode) -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumNonlinSolvIters" cARKodeGetNumNonlinSolvIters :: (SolverObject ARKode) -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumNonlinSolvConvFails" cARKodeGetNumNonlinSolvConvFails :: (SolverObject ARKode) -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumJacEvals" cARKodeGetNumJacEvals :: (SolverObject ARKode) -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKStepGetNumRhsEvals" cARKStepGetNumRhsEvals :: (SolverObject ARKode) -> Ptr CLong -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumLinRhsEvals" cARKodeGetNumLinRhsEvals :: (SolverObject ARKode) -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumStepAttempts" cARKodeGetNumStepAttempts :: (SolverObject ARKode) -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetNumGEvals" cARKodeGetNumGEvals :: (SolverObject ARKode) -> Ptr CLong -> IO (Flag ARKode)

foreign import ccall "ARKodeGetEstLocalErrors" cARKodeGetEstLocalErrors :: (SolverObject ARKode) -> N_Vector -> IO (Flag ARKode)

foreign import ccall "ARKodeGetErrWeights" cARKodeGetErrWeights :: (SolverObject ARKode) -> N_Vector -> IO (Flag ARKode)

data ARKode

instance IsLabel "SUCCESS" (Flag ARKode) where
  fromLabel = Flag ARK_SUCCESS

instance IsLabel "TOO_CLOSE" (Flag ARKode) where
  fromLabel = Flag ARK_TOO_CLOSE

instance IsLabel "ROOT_RETURN" (Flag ARKode) where
  fromLabel = Flag ARK_ROOT_RETURN
