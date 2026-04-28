{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Numeric.Sundials.Bindings.IDA where

import Foreign
import Foreign.C
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import GHC.OverloadedLabels

foreign import ccall "IDACreate" cIDACreate :: SUNContext -> IO (SolverObject IDA)

foreign import ccall "IDAFree" cIDAFree :: Ptr (SolverObject IDA) -> IO ()

foreign import ccall "IDAInit" cIDAInit :: (SolverObject IDA) -> FunPtr IDAResFn -> CDouble -> N_Vector -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDASetUserData" cIDASetUserData :: (SolverObject IDA) -> Ptr UserData -> IO (Flag IDA)

foreign import ccall "IDASetMinStep" cIDASetMinStep :: (SolverObject IDA) -> CDouble -> IO (Flag IDA)

foreign import ccall "IDASetMaxNumSteps" cIDASetMaxNumSteps :: (SolverObject IDA) -> SunIndexType -> IO (Flag IDA)

foreign import ccall "IDASetMaxErrTestFails" cIDASetMaxErrTestFails :: (SolverObject IDA) -> CInt -> IO (Flag IDA)

foreign import ccall "IDASVtolerances" cIDASVtolerances :: (SolverObject IDA) -> CDouble -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDARootInit" cIDARootInit :: (SolverObject IDA) -> CInt -> FunPtr IDARootFn -> IO (Flag IDA)

foreign import ccall "IDASetRootDirection" cIDASetRootDirection :: (SolverObject IDA) -> Ptr CInt -> IO (Flag IDA)

foreign import ccall "IDASetNoInactiveRootWarn" cIDASetNoInactiveRootWarn :: (SolverObject IDA) -> IO (Flag IDA)

foreign import ccall "IDASetLinearSolver" cIDASetLinearSolver :: (SolverObject IDA) -> SUNLinearSolver -> SUNMatrix -> IO (Flag IDA)

foreign import ccall "IDASetInitStep" cIDASetInitStep :: (SolverObject IDA) -> CDouble -> IO (Flag IDA)

foreign import ccall "IDASolve" cIDASolve :: (SolverObject IDA) -> CDouble -> Ptr CDouble -> N_Vector -> N_Vector -> CInt -> IO (Flag IDA)

foreign import ccall "IDAGetCurrentTime" cIDAGetCurrentTime :: (SolverObject IDA) -> Ptr CDouble -> IO (Flag IDA)

foreign import ccall "IDAReInit"
  cIDAReInit ::
    (SolverObject IDA) ->
    CDouble ->
    N_Vector ->
    N_Vector ->
    IO (Flag IDA)

foreign import ccall "IDAGetRootInfo" cIDAGetRootInfo :: (SolverObject IDA) -> Ptr CInt -> IO (Flag IDA)

foreign import ccall "IDASetJacFn" cIDASetJacFn :: (SolverObject IDA) -> FunPtr IDALsJacFn -> IO (Flag IDA)

foreign import ccall "IDAGetNumSteps" cIDAGetNumSteps :: (SolverObject IDA) -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumLinSolvSetups" cIDAGetNumLinSolvSetups :: (SolverObject IDA) -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumErrTestFails" cIDAGetNumErrTestFails :: (SolverObject IDA) -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumNonlinSolvIters" cIDAGetNumNonlinSolvIters :: (SolverObject IDA) -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumNonlinSolvConvFails" cIDAGetNumNonlinSolvConvFails :: (SolverObject IDA) -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumJacEvals" cIDAGetNumJacEvals :: (SolverObject IDA) -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDASetMaxStep" cIDASetMaxStep :: (SolverObject IDA) -> CDouble -> IO (Flag IDA)

foreign import ccall "IDAGetNumResEvals" cIDAGetNumResEvals :: (SolverObject IDA) -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumLinResEvals" cIDAGetNumLinResEvals :: (SolverObject IDA) -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetNumGEvals" cIDAGetNumGEvals :: (SolverObject IDA) -> Ptr CLong -> IO (Flag IDA)

foreign import ccall "IDAGetEstLocalErrors" cIDAGetEstLocalErrors :: (SolverObject IDA) -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDAGetErrWeights" cIDAGetErrWeights :: (SolverObject IDA) -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDACalcIC" cIDACalcIC :: (SolverObject IDA) -> CInt -> CDouble -> IO (Flag IDA)

foreign import ccall "IDAGetConsistentIC" cIDAGetConsistentIC :: (SolverObject IDA) -> N_Vector -> N_Vector -> IO (Flag IDA)

foreign import ccall "IDASetId" cIDASetId :: (SolverObject IDA) -> N_Vector -> IO (Flag IDA)

data IDA

instance IsLabel "SUCCESS" (Flag IDA) where
  fromLabel = Flag IDA_SUCCESS

instance IsLabel "TOO_CLOSE" (Flag IDA) where
  fromLabel = Flag IDA_TOO_CLOSE

instance IsLabel "ROOT_RETURN" (Flag IDA) where
  fromLabel = Flag IDA_ROOT_RETURN
