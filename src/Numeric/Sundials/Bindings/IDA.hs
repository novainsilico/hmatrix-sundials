{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Numeric.Sundials.Bindings.IDA where

import Data.Void
import Foreign
import Foreign.C
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign

-- | An opaque pointer to a IDAMem
newtype IDAMem = IDAMem (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "IDACreate" cIDACreate :: SUNContext -> IO IDAMem

foreign import ccall "IDAFree" cIDAFree :: Ptr IDAMem -> IO ()

foreign import ccall "IDAInit" cIDAInit :: IDAMem -> FunPtr IDAResFn -> CDouble -> N_Vector -> N_Vector -> IO CInt

foreign import ccall "IDASetUserData" cIDASetUserData :: IDAMem -> Ptr UserData -> IO CInt

foreign import ccall "IDASetMinStep" cIDASetMinStep :: IDAMem -> CDouble -> IO CInt

foreign import ccall "IDASetMaxNumSteps" cIDASetMaxNumSteps :: IDAMem -> SunIndexType -> IO CInt

foreign import ccall "IDASetMaxErrTestFails" cIDASetMaxErrTestFails :: IDAMem -> CInt -> IO CInt

foreign import ccall "IDASVtolerances" cIDASVtolerances :: IDAMem -> CDouble -> N_Vector -> IO CInt

foreign import ccall "IDARootInit" cIDARootInit :: IDAMem -> CInt -> FunPtr IDARootFn -> IO CInt

foreign import ccall "IDASetRootDirection" cIDASetRootDirection :: IDAMem -> Ptr CInt -> IO CInt

foreign import ccall "IDASetNoInactiveRootWarn" cIDASetNoInactiveRootWarn :: IDAMem -> IO CInt

foreign import ccall "IDASetLinearSolver" cIDASetLinearSolver :: IDAMem -> SUNLinearSolver -> SUNMatrix -> IO CInt

foreign import ccall "IDASetInitStep" cIDASetInitStep :: IDAMem -> CDouble -> IO CInt

foreign import ccall "IDASolve" cIDASolve :: IDAMem -> CDouble -> Ptr CDouble -> N_Vector -> N_Vector -> CInt -> IO CInt

foreign import ccall "IDAGetCurrentTime" cIDAGetCurrentTime :: IDAMem -> Ptr CDouble -> IO CInt

foreign import ccall "IDAReInit"
  cIDAReInit ::
    IDAMem ->
    CDouble ->
    N_Vector ->
    N_Vector ->
    IO CInt

foreign import ccall "IDAGetRootInfo" cIDAGetRootInfo :: IDAMem -> Ptr CInt -> IO CInt

foreign import ccall "IDASetJacFn" cIDASetJacFn :: IDAMem -> FunPtr IDALsJacFn -> IO CInt

foreign import ccall "IDAGetNumSteps" cIDAGetNumSteps :: IDAMem -> Ptr CLong -> IO CInt

foreign import ccall "IDAGetNumLinSolvSetups" cIDAGetNumLinSolvSetups :: IDAMem -> Ptr CLong -> IO CInt

foreign import ccall "IDAGetNumErrTestFails" cIDAGetNumErrTestFails :: IDAMem -> Ptr CLong -> IO CInt

foreign import ccall "IDAGetNumNonlinSolvIters" cIDAGetNumNonlinSolvIters :: IDAMem -> Ptr CLong -> IO CInt

foreign import ccall "IDAGetNumNonlinSolvConvFails" cIDAGetNumNonlinSolvConvFails :: IDAMem -> Ptr CLong -> IO CInt

foreign import ccall "IDAGetNumJacEvals" cIDAGetNumJacEvals :: IDAMem -> Ptr CLong -> IO CInt

foreign import ccall "IDASetMaxStep" cIDASetMaxStep :: IDAMem -> CDouble -> IO CInt

foreign import ccall "IDAGetNumResEvals" cIDAGetNumResEvals :: IDAMem -> Ptr CLong -> IO CInt

foreign import ccall "IDAGetNumLinResEvals" cIDAGetNumLinResEvals :: IDAMem -> Ptr CLong -> IO CInt

foreign import ccall "IDAGetNumGEvals" cIDAGetNumGEvals :: IDAMem -> Ptr CLong -> IO CInt

foreign import ccall "IDAGetEstLocalErrors" cIDAGetEstLocalErrors :: IDAMem -> N_Vector -> IO CInt

foreign import ccall "IDAGetErrWeights" cIDAGetErrWeights :: IDAMem -> N_Vector -> IO CInt

foreign import ccall "IDACalcIC" cIDACalcIC :: IDAMem -> CInt -> CDouble -> IO CInt

foreign import ccall "IDAGetConsistentIC" cIDAGetConsistentIC :: IDAMem -> N_Vector -> N_Vector -> IO CInt

foreign import ccall "IDASetId" cIDASetId :: IDAMem -> N_Vector -> IO CInt
