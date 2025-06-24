{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

-- | This module contains most bindings relevant to sundial and utilities
-- (matrix, vector) which are used by all solvers.
module Numeric.Sundials.Bindings.Sundials where

import Control.Exception
import Control.Monad (when)
import Data.Void
import Foreign
import Foreign.C
import GHC.Stack
import Numeric.Sundials.Common
import Numeric.Sundials.Foreign

-- * Sundial context

newtype SUNContext = SUNContext (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "SUNContext_Create" cSUNContext_Create :: Int -> Ptr SUNContext -> IO CInt

{-# WARNING cSUNContext_Create "Prefere 'withSUNContext' to avoid memory leaks." #-}

foreign import ccall "SUNContext_Free" cSUNContext_Free :: Ptr SUNContext -> IO CInt

withSUNContext :: (HasCallStack) => (SUNContext -> IO a) -> IO a
withSUNContext cont = do
  alloca $ \ptr -> do
    let create = do
          cSUNContext_Create 0 ptr >>= failOnError
        destroy = cSUNContext_Free ptr
    bracket_ create destroy $ do
      sunctx <- peek ptr
      cont sunctx

-- * Error handling

foreign import ccall "SUNGetErrMsg" cSUNGetErrMsg :: CInt -> IO CString

-- | Push an error to the stack
--
-- Be careful, you may need to first call cSUNContext_ClearErrHandlers to
-- remove all previous handler
foreign import ccall "SUNContext_PushErrHandler" cSUNContext_PushErrHandler :: SUNContext -> FunPtr ReportErrorFnNew -> Ptr () -> IO CInt

{-# WARNING cSUNContext_PushErrHandler "This function is dangerous because it adds an handler on top of previously existing handlers. It can hence results in duplicated logging messages as well as extra (e.g. not filtered by the specified handler). Be sure to call cSUNContext_ClearErrHandlers' before. Or call setErrorHandler" #-}

foreign import ccall "SUNContext_ClearErrHandlers" cSUNContext_ClearErrHandlers :: SUNContext -> IO CInt

failOnError :: CInt -> IO ()
failOnError errCode = do
  when (errCode /= 0) $ do
    errMsg <- cSUNGetErrMsg errCode
    msg <- peekCString errMsg
    error msg


-- | Set an error handler (and dispable the other ones)
setErrorHandler :: SUNContext -> FunPtr ReportErrorFnNew -> IO ()
setErrorHandler sunctx handler = do
  cSUNContext_ClearErrHandlers sunctx >>= failOnError
  cSUNContext_PushErrHandler sunctx handler nullPtr >>= failOnError

-- * Vectors

-- | An opaque pointer to an N_Vector
newtype N_Vector = N_Vector {getSunVector :: Ptr SunVector}

foreign import ccall "N_VNew_Serial" cN_VNew_Serial :: SunIndexType -> SUNContext -> IO N_Vector

{-# WARNING cN_VNew_Serial "Prefere 'withNVector_Serial' to avoid memory leaks." #-}

foreign import ccall "N_VDestroy" cN_VDestroy :: N_Vector -> IO ()


-- TODO: build a wrapper which init the vector
withNVector_Serial :: SunIndexType -> SUNContext -> Int -> (N_Vector -> IO a) -> IO a
withNVector_Serial t suncontext errCode f = do
  let create = do
        res@(N_Vector ptr) <- cN_VNew_Serial t suncontext
        when (ptr == nullPtr) $ do
          throwIO $ ReturnCodeWithMessage "Failure in N_VNew_Serial" errCode
        pure res
  bracket create cN_VDestroy f

-- | Set element of a vector
cNV_Ith_S :: N_Vector -> Int -> CDouble -> IO ()
cNV_Ith_S (N_Vector ptr) i v = do
  qtr <- getContentPtr ptr
  rtr <- getData qtr
  pokeElemOff rtr i v

-- | Get element of a vector
cNV_Ith_S' :: N_Vector -> Int -> IO CDouble
cNV_Ith_S' (N_Vector ptr) i = do
  qtr <- getContentPtr ptr
  rtr <- getData qtr
  peekElemOff rtr i

foreign import ccall "N_VGetArrayPointer" cN_VGetArrayPointer :: N_Vector -> IO (Ptr CDouble)

-- * Matrix

-- | Opaque structure for a matrix
newtype SUNMatrix = SUNMatrix (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "SUNSparseMatrix" cSUNSparseMatrix :: SunIndexType -> SunIndexType -> CInt -> CInt -> SUNContext -> IO SUNMatrix

{-# WARNING cSUNSparseMatrix "Prefere 'withSUNSparseMatrix' to avoid memory leaks" #-}

foreign import ccall "SUNDenseMatrix" cSUNDenseMatrix :: SunIndexType -> SunIndexType -> SUNContext -> IO SUNMatrix

{-# WARNING cSUNDenseMatrix "Prefere 'withSUNDenseMatrix' to avoid memory leaks" #-}

foreign import ccall "SUNMatDestroy" cSUNMatDestroy :: SUNMatrix -> IO ()

withSUNDenseMatrix :: (HasCallStack) => SunIndexType -> SunIndexType -> SUNContext -> CInt -> (SUNMatrix -> IO a) -> IO a
withSUNDenseMatrix dim dim' sunctx errCode f = do
  let create = do
        mat@(SUNMatrix ptr) <- cSUNDenseMatrix dim dim' sunctx

        when (ptr == nullPtr) $ do
          throwIO $ ReturnCode $ fromIntegral errCode

        pure mat
  bracket create cSUNMatDestroy f

withSUNSparseMatrix :: (HasCallStack) => SunIndexType -> SunIndexType -> CInt -> CInt -> SUNContext -> CInt -> (SUNMatrix -> IO a) -> IO a
withSUNSparseMatrix dim dim' jac set sunctx errCode f = do
  let create = do
        mat@(SUNMatrix ptr) <- cSUNSparseMatrix dim dim' jac set sunctx

        when (ptr == nullPtr) $ do
          throwIO $ ReturnCode $ fromIntegral errCode

        pure mat
  bracket create cSUNMatDestroy f

-- * Solvers

newtype SUNLinearSolver = SUNLinearSolver (Ptr Void)
  deriving newtype (Storable)

foreign import ccall "SUNLinSol_KLU" cSUNLinSol_KLU :: N_Vector -> SUNMatrix -> SUNContext -> IO SUNLinearSolver

{-# WARNING cSUNLinSol_KLU "Prefer 'withSUNLinSol_KLU' to avoid memory leak." #-}

foreign import ccall "SUNLinSol_Dense" cSUNLinSol_Dense :: N_Vector -> SUNMatrix -> SUNContext -> IO SUNLinearSolver

{-# WARNING cSUNLinSol_Dense "Prefer 'withSUNLinSol_Dense' to avoid memory leak." #-}

foreign import ccall "SUNLinSolFree" cSUNLinSolFree :: SUNLinearSolver -> IO CInt

withSUNLinSol_Dense :: (HasCallStack) => N_Vector -> SUNMatrix -> SUNContext -> CInt -> (SUNLinearSolver -> IO a) -> IO a
withSUNLinSol_Dense vec mat sunctx errCode f = do
  let create = do
        ls@(SUNLinearSolver ptr) <- cSUNLinSol_Dense vec mat sunctx

        when (ptr == nullPtr) $ do
          throwIO $ ReturnCode $ fromIntegral errCode

        pure ls
  bracket create cSUNLinSolFree f

withSUNLinSol_KLU :: (HasCallStack) => N_Vector -> SUNMatrix -> SUNContext -> CInt -> (SUNLinearSolver -> IO a) -> IO a
withSUNLinSol_KLU vec mat sunctx errCode f = do
  let create = do
        ls@(SUNLinearSolver ptr) <- cSUNLinSol_KLU vec mat sunctx

        when (ptr == nullPtr) $ do
          throwIO $ ReturnCode $ fromIntegral errCode

        pure ls
  bracket create cSUNLinSolFree f

-- * Solver loop logic

-- Tries to copy the previous semantic in C when we had:
data ReturnCode
  = -- return, with error code, skipping finish (diagnostics and cleanup)
    ReturnCode Int
  | -- return with log message, skipping finish (diags and cleanup)
    ReturnCodeWithMessage String Int
  | -- "goto" finish (do diag and cleanup)
    Finish LoopState
  | -- break the loop, hence go to finish (do diag and cleanup)
    Break LoopState
  deriving (Exception, Show)

-- | The loop state, most could be just carried by recursive function call, but
-- the C semantic was kinda interleaved, so doing that for now.
data LoopState = LoopState
  { -- output_ind tracks the current row into the c_output_mat matrix.
    -- if differs from input_ind because of the extra rows corresponding to events.
    output_ind :: {-# UNPACK #-} !Int,
    -- input_ind tracks the current index into the c_sol_time array
    input_ind :: {-# UNPACK #-} !Int,
    -- event_ind tracks the current event number
    event_ind :: {-# UNPACK #-} !Int,
    t_start :: {-# UNPACK #-} !CDouble
  }
  deriving (Show)

{- NOTE [SAFETY]
   ~~~~~~~~~~~~~~~~~~~

All the call are marked (implicitly) as "safe" because they can, in theory, all call the log callback and hence "unsafe" is not relevant.

-}
