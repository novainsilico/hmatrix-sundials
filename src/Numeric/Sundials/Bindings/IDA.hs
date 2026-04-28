{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Numeric.Sundials.Bindings.IDA where

import Data.Void
import Foreign
import Foreign.C
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import GHC.OverloadedLabels

-- | An opaque pointer to a IDAMem
newtype IDAMem = IDAMem (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "IDACreate" cIDACreate :: SUNContext -> IO IDAMem

foreign import ccall "IDAFree" cIDAFree :: Ptr IDAMem -> IO ()

foreign import ccall "IDAInit" cIDAInit :: IDAMem -> FunPtr IDAResFn -> CDouble -> N_Vector -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDASetUserData" cIDASetUserData :: IDAMem -> Ptr UserData -> IO (Flag IDA)

foreign import ccall "IDASetMinStep" cIDASetMinStep :: IDAMem -> CDouble -> IO (Flag IDA)

foreign import ccall "IDASetMaxNumSteps" cIDASetMaxNumSteps :: IDAMem -> SunIndexType -> IO (Flag IDA)

foreign import ccall "IDASetMaxErrTestFails" cIDASetMaxErrTestFails :: IDAMem -> CInt -> IO (Flag IDA)

foreign import ccall "IDASVtolerances" cIDASVtolerances :: IDAMem -> CDouble -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDARootInit" cIDARootInit :: IDAMem -> CInt -> FunPtr IDARootFn -> IO (Flag IDA)

foreign import ccall "IDASetRootDirection" cIDASetRootDirection :: IDAMem -> Ptr CInt -> IO (Flag IDA)

foreign import ccall "IDASetNoInactiveRootWarn" cIDASetNoInactiveRootWarn :: IDAMem -> IO (Flag IDA)

foreign import ccall "IDASetLinearSolver" cIDASetLinearSolver :: IDAMem -> SUNLinearSolver -> SUNMatrix -> IO (Flag IDA)

foreign import ccall "IDASetInitStep" cIDASetInitStep :: IDAMem -> CDouble -> IO (Flag IDA)

foreign import ccall "IDASolve" cIDASolve :: IDAMem -> CDouble -> Ptr CDouble -> N_Vector -> N_Vector -> CInt -> IO (Flag IDA)

foreign import ccall "IDAGetCurrentTime" cIDAGetCurrentTime :: IDAMem -> Ptr CDouble -> IO (Flag IDA)

foreign import ccall "IDAReInit"
  cIDAReInit ::
    IDAMem ->
    CDouble ->
    N_Vector ->
    N_Vector ->
    IO (Flag IDA)

foreign import ccall "IDAGetRootInfo" cIDAGetRootInfo :: IDAMem -> Ptr CInt -> IO (Flag IDA)

foreign import ccall "IDASetJacFn" cIDASetJacFn :: IDAMem -> FunPtr IDALsJacFn -> IO (Flag IDA)

foreign import ccall "IDAGetNumSteps" cIDAGetNumSteps :: IDAMem -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumLinSolvSetups" cIDAGetNumLinSolvSetups :: IDAMem -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumErrTestFails" cIDAGetNumErrTestFails :: IDAMem -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumNonlinSolvIters" cIDAGetNumNonlinSolvIters :: IDAMem -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumNonlinSolvConvFails" cIDAGetNumNonlinSolvConvFails :: IDAMem -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumJacEvals" cIDAGetNumJacEvals :: IDAMem -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDASetMaxStep" cIDASetMaxStep :: IDAMem -> CDouble -> IO (Flag IDA)

foreign import ccall "IDAGetNumResEvals" cIDAGetNumResEvals :: IDAMem -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumLinResEvals" cIDAGetNumLinResEvals :: IDAMem -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumGEvals" cIDAGetNumGEvals :: IDAMem -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetEstLocalErrors" cIDAGetEstLocalErrors :: IDAMem -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDAGetErrWeights" cIDAGetErrWeights :: IDAMem -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDACalcIC" cIDACalcIC :: IDAMem -> CInt -> CDouble -> IO (Flag IDA)

foreign import ccall "IDAGetConsistentIC" cIDAGetConsistentIC :: IDAMem -> N_Vector -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDASetId" cIDASetId :: IDAMem -> N_Vector -> IO (Flag IDA)

data IDA

instance IsLabel "SUCCESS" (Flag IDA) where
  fromLabel = Flag IDA_SUCCESS

instance IsLabel "TOO_CLOSE" (Flag IDA) where
  fromLabel = Flag IDA_TOO_CLOSE

instance IsLabel "ROOT_RETURN" (Flag IDA) where
  fromLabel = Flag IDA_ROOT_RETURN
