{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Numeric.Sundials.Bindings.CVode where

import Foreign
import Foreign.C
import Numeric.Sundials.Bindings.Sundials
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign
import GHC.OverloadedLabels

foreign import ccall "CVodeCreate" cCVodeCreate :: CInt -> SUNContext -> IO (SolverObject CVode)

foreign import ccall "CVodeFree" cCVodeFree :: Ptr (SolverObject CVode) -> IO ()

foreign import ccall "CVodeInit" cCVodeInit :: (SolverObject CVode) -> FunPtr OdeRhsCType -> SunRealType -> N_Vector -> IO (Flag CVode)

foreign import ccall "CVodeSetUserData" cCVodeSetUserData :: (SolverObject CVode) -> Ptr UserData -> IO (Flag CVode)

foreign import ccall "CVodeSetMinStep" cCVodeSetMinStep :: (SolverObject CVode) -> CDouble -> IO (Flag CVode)

foreign import ccall "CVodeSetMaxStep" cCVodeSetMaxStep :: (SolverObject CVode) -> CDouble -> IO (Flag CVode)

foreign import ccall "CVodeSetMaxNumSteps" cCVodeSetMaxNumSteps :: (SolverObject CVode) -> SunIndexType -> IO (Flag CVode)

foreign import ccall "CVodeSetMaxErrTestFails" cCVodeSetMaxErrTestFails :: (SolverObject CVode) -> CInt -> IO (Flag CVode)

foreign import ccall "CVodeSVtolerances" cCVodeSVtolerances :: (SolverObject CVode) -> CDouble -> N_Vector -> IO (Flag CVode)

foreign import ccall "CVodeRootInit" cCVodeRootInit :: (SolverObject CVode) -> CInt -> FunPtr EventConditionCType -> IO (Flag CVode)

foreign import ccall "CVodeSetRootDirection" cCVodeSetRootDirection :: (SolverObject CVode) -> Ptr CInt -> IO (Flag CVode)

foreign import ccall "CVodeSetNoInactiveRootWarn" cCVodeSetNoInactiveRootWarn :: (SolverObject CVode) -> IO (Flag CVode)

foreign import ccall "CVodeSetLinearSolver" cCVodeSetLinearSolver :: (SolverObject CVode) -> SUNLinearSolver -> SUNMatrix -> IO (Flag CVode)

foreign import ccall "CVodeSetInitStep" cCVodeSetInitStep :: (SolverObject CVode) -> CDouble -> IO (Flag CVode)

foreign import ccall "CVode" cCVode :: (SolverObject CVode) -> CDouble -> N_Vector -> Ptr CDouble -> CInt -> IO (Flag CVode)

foreign import ccall "CVodeReInit" cCVodeReInit :: (SolverObject CVode) -> CDouble -> N_Vector -> IO ()

foreign import ccall "CVodeGetRootInfo" cCVodeGetRootInfo :: (SolverObject CVode) -> Ptr CInt -> IO (Flag CVode)

foreign import ccall "CVodeSetJacFn" cCVodeSetJacFn :: (SolverObject CVode) -> FunPtr OdeJacobianCType -> IO (Flag CVode)

foreign import ccall "CVodeGetNumSteps" cCVodeGetNumSteps :: (SolverObject CVode) -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumLinSolvSetups" cCVodeGetNumLinSolvSetups :: (SolverObject CVode) -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumErrTestFails" cCVodeGetNumErrTestFails :: (SolverObject CVode) -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumNonlinSolvIters" cCVodeGetNumNonlinSolvIters :: (SolverObject CVode) -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumNonlinSolvConvFails" cCVodeGetNumNonlinSolvConvFails :: (SolverObject CVode) -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumJacEvals" cCVodeGetNumJacEvals :: (SolverObject CVode) -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumRhsEvals" cCVodeGetNumRhsEvals :: (SolverObject CVode) -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumLinRhsEvals" cCVodeGetNumLinRhsEvals :: (SolverObject CVode) -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetNumGEvals" cCVodeGetNumGEvals :: (SolverObject CVode) -> Ptr CLong -> IO (Flag CVode)

foreign import ccall "CVodeGetEstLocalErrors" cCVodeGetEstLocalErrors :: (SolverObject CVode) -> N_Vector -> IO (Flag CVode)

foreign import ccall "CVodeGetErrWeights" cCVodeGetErrWeights :: (SolverObject CVode) -> N_Vector -> IO (Flag CVode)

data CVode

instance IsLabel "SUCCESS" (Flag CVode) where
  fromLabel = Flag CV_SUCCESS

instance IsLabel "TOO_CLOSE" (Flag CVode) where
  fromLabel = Flag CV_TOO_CLOSE

instance IsLabel "ROOT_RETURN" (Flag CVode) where
  fromLabel = Flag CV_ROOT_RETURN
