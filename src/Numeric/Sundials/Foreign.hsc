{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE EmptyDataDecls #-}
{-# LANGUAGE PatternSynonyms #-}

module Numeric.Sundials.Foreign
  ( getDataFromContents
  , putDataInContents
  , cV_ADAMS
  , cV_BDF
  , vectorToC
  , SunIndexType
  , SunRealType
  , SunMatrix(..)
  , SparsePattern(..)
  , SparseMatrix(..)
  , SunVector(..)
    -- * Offsets
    -- ** NVector
  , nvectorContentOffset
    -- ** NVector_SERIAL
  , nvectorContentSerialLengthOffset
  , nvectorContentSerialDataOffset
    -- ** SUNMatrix
  , sunmatrixContentOffset
    -- ** SUNMatrix_DENSE
  , sunmatrixContentDenseDataOffset
    -- ** SUNMatrix_SPARSE
  , sunmatrixContentSparseIndexvalsOffset
  , sunmatrixContentSparseIndexptrsOffset
  , sunmatrixContentSparseDataOffset
  , sunmatrixContentSparseNnzOffset
    -- * Methods
  , bACKWARD_EULER_1_1
  , iMPLICIT_MIDPOINT_1_2
  , aRK2_DIRK_3_1_2
  , eSDIRK325L2SA_5_2_3
  , eSDIRK436L2SA_6_3_4
  , fORWARD_EULER_1_1
  , bOGACKI_SHAMPINE_4_2_3
  , sOFRONIOU_SPALETTA_5_3_4
  , tSITOURAS_7_4_5
  , vERNER_10_6_7
  , vERNER_16_8_9
  , pattern CSC_MAT

  , pattern CV_NORMAL
  , pattern CV_ROOT_RETURN
  , pattern CV_SUCCESS
  , pattern CV_TOO_CLOSE

  , pattern ARK_NORMAL
  , pattern ARK_ROOT_RETURN
  , pattern ARK_SUCCESS
  , pattern ARK_TOO_CLOSE

  , pattern IDA_NORMAL
  , pattern IDA_ROOT_RETURN
  , pattern IDA_SUCCESS
  , pattern IDA_ILL_INPUT
  , pattern IDA_Y_INIT
  , pattern IDA_YA_YDP_INIT
  , pattern IDA_ERR_FAIL
  , pattern IDA_CONV_FAIL
  , pattern IDA_TOO_MUCH_WORK
  , pattern IDA_TOO_CLOSE

  , getContentPtr, getData
  ) where

import           Foreign
import           Foreign.C.Types

import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VM
import Data.IORef
import Data.Coerce
import Control.Monad
import Control.Exception
import Text.Printf (printf)

#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <arkode/arkode.h>
#include <arkode/arkode_arkstep.h>
#include <cvode/cvode.h>
#include <ida/ida.h>

-- A version of #type that produces CDouble etc.
-- https://github.com/haskell/hsc2hs/issues/51
#define hsc_ctype(t...)                                     \
    if ((t)(int)(t)1.4 == (t)1.4)                           \
        hsc_printf ("%s%lu",                                \
                (t)(-1) < (t)0 ? "Int" : "Word",            \
                (unsigned long)sizeof (t) * 8);             \
    else                                                    \
        hsc_printf ("%s",                                   \
                sizeof (t) >  sizeof (double) ? "LDouble" : \
                sizeof (t) == sizeof (double) ? "CDouble" : \
                "CFloat");

data SunVector = SunVector
  { sunVecN    :: SunIndexType
  , sunVecVals :: VS.Vector CDouble
  }

data SunMatrix = SunMatrix
  { rows :: SunIndexType
  , cols :: SunIndexType
  , vals :: VS.Vector CDouble
    -- ^ matrix entries in the column-major order
  }

-- | A sparse pattern: a column-wise stored matrix that has 0 for zero
-- entries and 1 for non-zero entries.
newtype SparsePattern = SparsePattern (VS.Vector Int8)
  deriving Show

-- | A sparse matrix
data SparseMatrix = SparseMatrix SparsePattern SunMatrix

type SunIndexType = #ctype sunindextype
type SunRealType = #ctype sunrealtype

getMatrixDataFromContents :: Ptr SunMatrix -> IO SunMatrix
getMatrixDataFromContents ptr = do
  qtr <- getContentMatrixPtr ptr
  rs  <- getNRows qtr
  cs  <- getNCols qtr
  rtr <- getMatrixData qtr
  vs  <- vectorFromC (fromIntegral $ rs * cs) rtr
  return $ SunMatrix { rows = rs, cols = cs, vals = vs }

putMatrixDataFromContents :: SunMatrix -> Ptr SunMatrix -> IO ()
putMatrixDataFromContents mat ptr = do
  let rs = rows mat
      cs = cols mat
      vs = vals mat
  qtr <- getContentMatrixPtr ptr
  putNRows rs qtr
  putNCols cs qtr
  rtr <- getMatrixData qtr
  vectorToC vs (fromIntegral $ rs * cs) rtr

instance Storable SunVector where
  poke p v    = putDataInContents (sunVecVals v) (fromIntegral $ sunVecN v) p
  peek p      = do (l, v) <- getDataFromContents p
                   return $ SunVector { sunVecN = fromIntegral l
                                      , sunVecVals = v
                                      }
  sizeOf _    = error "sizeOf not supported for SunVector"
  alignment _ = error "alignment not supported for SunVector"

instance Storable SunMatrix where
  poke        = flip putMatrixDataFromContents
  peek        = getMatrixDataFromContents
  sizeOf _    = error "sizeOf not supported for SunMatrix"
  alignment _ = error "alignment not supported for SunMatrix"

instance Storable SparseMatrix where
  poke p (SparseMatrix (SparsePattern spat) SunMatrix{..}) = do
    content_ptr <- getContentMatrixPtr p
    sparse_type :: CInt <- #{peek struct _SUNMatrixContent_Sparse, sparsetype} content_ptr
    unless (sparse_type == cSC_MAT) $
      throwIO . ErrorCall $ printf "Sparse SUNMatrix poke: only CSC matrices are supported, got %d" (fromIntegral sparse_type :: Int)
    #{poke struct _SUNMatrixContent_Sparse, M} content_ptr rows
    #{poke struct _SUNMatrixContent_Sparse, N} content_ptr cols
    indexvals :: Ptr SunIndexType <- #{peek struct _SUNMatrixContent_Sparse, indexvals} content_ptr
    indexptrs :: Ptr SunIndexType <- #{peek struct _SUNMatrixContent_Sparse, indexptrs} content_ptr
    data_ :: Ptr SunRealType <- #{peek struct _SUNMatrixContent_Sparse, data} content_ptr
    max_entries :: SunIndexType <- #{peek struct _SUNMatrixContent_Sparse, NNZ} content_ptr

    when (VS.any (`notElem` [0,1]) spat) $
      throwIO $ ErrorCall $ "Illegal sparse pattern: " ++ show spat
    when (max_entries < VS.sum (VS.map fromIntegral spat)) $
      throwIO $ ErrorCall $ "Not enough space in sparse matrix for the sparse pattern"

    out_ix_ref <- newIORef 0
    forM_ [0 .. cols - 1] $ \cur_col -> do
      readIORef out_ix_ref >>= \out_ix ->
        pokeElemOff indexptrs (fromIntegral cur_col) (fromIntegral out_ix)
      forM_ [0 .. rows - 1] $ \cur_row -> do
        let
          in_ix = fromIntegral $ cur_row + cur_col * rows
          non_zero = (spat VS.! in_ix) /= 0
          val = vals VS.! in_ix
        -- Make sure the diagonal entries are always allocated, to be able
        -- to store the I - γJ matrix.
        when (non_zero || cur_row == cur_col) $ do
          out_ix <- readIORef out_ix_ref
          writeIORef out_ix_ref $! out_ix+1
          pokeElemOff data_ out_ix (coerce val)
          pokeElemOff indexvals out_ix (fromIntegral cur_row)
    readIORef out_ix_ref >>= \out_ix ->
      pokeElemOff indexptrs (fromIntegral cols) (fromIntegral out_ix)


  peek        = error "peek not supported for a sparse SunMatrix"
  sizeOf _    = error "sizeOf not supported for SunMatrix"
  alignment _ = error "alignment not supported for SunMatrix"

vectorFromC :: Storable a => Int -> Ptr a -> IO (VS.Vector a)
vectorFromC len ptr = do
  ptr' <- newForeignPtr_ ptr
  VS.freeze $ VM.unsafeFromForeignPtr0 ptr' len

vectorToC :: Storable a => VS.Vector a -> Int -> Ptr a -> IO ()
vectorToC vec len ptr = do
  ptr' <- newForeignPtr_ ptr
  VS.copy (VM.unsafeFromForeignPtr0 ptr' len) vec

getDataFromContents :: Ptr SunVector -> IO (SunIndexType, VS.Vector CDouble)
getDataFromContents ptr = do
  qtr <- getContentPtr ptr
  rtr <- getData qtr
  len' <- getLength qtr
  v <- vectorFromC (fromIntegral len') rtr
  return (len', v)

putDataInContents :: VS.Vector CDouble -> Int -> Ptr SunVector -> IO ()
putDataInContents vec len ptr = do
  qtr <- getContentPtr ptr
  rtr <- getData qtr
  putLength (fromIntegral len) qtr
  vectorToC vec len rtr

----------------------------------------------------------------------
--                           Offsets
----------------------------------------------------------------------

nvectorContentOffset :: Num a => a
nvectorContentOffset = #offset struct _generic_N_Vector, content

nvectorContentSerialLengthOffset :: Num a => a
nvectorContentSerialLengthOffset = #offset struct _N_VectorContent_Serial, length

nvectorContentSerialDataOffset :: Num a => a
nvectorContentSerialDataOffset = #offset struct _N_VectorContent_Serial, data

sunmatrixContentOffset :: Num a => a
sunmatrixContentOffset = #offset struct _generic_SUNMatrix, content

sunmatrixContentDenseDataOffset :: Num a => a
sunmatrixContentDenseDataOffset = #offset struct _SUNMatrixContent_Dense, data

sunmatrixContentSparseIndexvalsOffset :: Num a => a
sunmatrixContentSparseIndexvalsOffset = #offset struct _SUNMatrixContent_Sparse, indexvals

sunmatrixContentSparseIndexptrsOffset :: Num a => a
sunmatrixContentSparseIndexptrsOffset = #offset struct _SUNMatrixContent_Sparse, indexptrs

sunmatrixContentSparseDataOffset :: Num a => a
sunmatrixContentSparseDataOffset = #offset struct _SUNMatrixContent_Sparse, data

sunmatrixContentSparseNnzOffset :: Num a => a
sunmatrixContentSparseNnzOffset = #offset struct _SUNMatrixContent_Sparse, NNZ



getContentMatrixPtr :: Storable a => Ptr b -> IO a
getContentMatrixPtr ptr = (#peek struct _generic_SUNMatrix, content) ptr

getNRows :: Ptr b -> IO SunIndexType
getNRows ptr = (#peek struct _SUNMatrixContent_Dense, M) ptr
putNRows :: SunIndexType -> Ptr b -> IO ()
putNRows nr ptr = (#poke struct _SUNMatrixContent_Dense, M) ptr nr

getNCols :: Ptr b -> IO SunIndexType
getNCols ptr = (#peek struct _SUNMatrixContent_Dense, N) ptr
putNCols :: SunIndexType -> Ptr b -> IO ()
putNCols nc ptr = (#poke struct _SUNMatrixContent_Dense, N) ptr nc

getMatrixData :: Storable a => Ptr b -> IO a
getMatrixData ptr = (#peek struct _SUNMatrixContent_Dense, data) ptr

getContentPtr :: Storable a => Ptr b -> IO a
getContentPtr ptr = (#peek struct _generic_N_Vector, content) ptr

getData :: Storable a => Ptr b -> IO a
getData ptr = (#peek struct _N_VectorContent_Serial, data) ptr

getLength :: Ptr b -> IO SunIndexType
getLength ptr = (#peek struct _N_VectorContent_Serial, length) ptr

putLength :: SunIndexType -> Ptr b -> IO ()
putLength l ptr = (#poke struct _N_VectorContent_Serial, length) ptr l

cSC_MAT :: CInt
cSC_MAT = #const CSC_MAT

cV_ADAMS :: CInt
cV_ADAMS = #const CV_ADAMS
cV_BDF :: CInt
cV_BDF = #const CV_BDF

-- Butcher table accessors
bACKWARD_EULER_1_1 :: CInt
bACKWARD_EULER_1_1 = #const ARKODE_BACKWARD_EULER_1_1
iMPLICIT_MIDPOINT_1_2 :: CInt
iMPLICIT_MIDPOINT_1_2 = #const ARKODE_IMPLICIT_MIDPOINT_1_2
aRK2_DIRK_3_1_2 :: CInt
aRK2_DIRK_3_1_2 = #const ARKODE_ARK2_DIRK_3_1_2
eSDIRK325L2SA_5_2_3 :: CInt
eSDIRK325L2SA_5_2_3 = #const ARKODE_ESDIRK325L2SA_5_2_3
eSDIRK436L2SA_6_3_4 :: CInt
eSDIRK436L2SA_6_3_4 = #const ARKODE_ESDIRK436L2SA_6_3_4
fORWARD_EULER_1_1 :: CInt
fORWARD_EULER_1_1 = #const ARKODE_FORWARD_EULER_1_1
bOGACKI_SHAMPINE_4_2_3 :: CInt
bOGACKI_SHAMPINE_4_2_3 = #const ARKODE_BOGACKI_SHAMPINE_4_2_3
sOFRONIOU_SPALETTA_5_3_4 :: CInt
sOFRONIOU_SPALETTA_5_3_4 = #const ARKODE_SOFRONIOU_SPALETTA_5_3_4
tSITOURAS_7_4_5 :: CInt
tSITOURAS_7_4_5 = #const ARKODE_TSITOURAS_7_4_5
vERNER_10_6_7 :: CInt
vERNER_10_6_7 = #const ARKODE_VERNER_10_6_7
vERNER_16_8_9 :: CInt
vERNER_16_8_9 = #const ARKODE_VERNER_16_8_9


pattern CSC_MAT :: CInt
pattern CSC_MAT = #const CSC_MAT

pattern CV_NORMAL, CV_SUCCESS, CV_ROOT_RETURN, CV_TOO_CLOSE :: CInt
pattern CV_NORMAL = #const CV_NORMAL
pattern CV_SUCCESS = #const CV_SUCCESS
pattern CV_ROOT_RETURN = #const CV_ROOT_RETURN
pattern CV_TOO_CLOSE = #const CV_TOO_CLOSE

pattern ARK_NORMAL, ARK_SUCCESS, ARK_ROOT_RETURN, ARK_TOO_CLOSE :: CInt
pattern ARK_NORMAL = #const ARK_NORMAL
pattern ARK_SUCCESS = #const ARK_SUCCESS
pattern ARK_ROOT_RETURN = #const ARK_ROOT_RETURN
pattern ARK_TOO_CLOSE = #const ARK_TOO_CLOSE

pattern IDA_NORMAL, IDA_SUCCESS, IDA_ROOT_RETURN, IDA_ILL_INPUT, IDA_Y_INIT, IDA_YA_YDP_INIT, IDA_TOO_CLOSE :: CInt
pattern IDA_NORMAL = #const IDA_NORMAL
pattern IDA_SUCCESS = #const IDA_SUCCESS
pattern IDA_ROOT_RETURN = #const IDA_ROOT_RETURN
pattern IDA_ILL_INPUT = #const IDA_ILL_INPUT
pattern IDA_Y_INIT = #const IDA_Y_INIT
pattern IDA_YA_YDP_INIT = #const IDA_YA_YDP_INIT
pattern IDA_ERR_FAIL = #const IDA_ERR_FAIL
pattern IDA_CONV_FAIL = #const IDA_CONV_FAIL
pattern IDA_TOO_MUCH_WORK = #const IDA_TOO_MUCH_WORK
pattern IDA_TOO_CLOSE = #const IDA_TOO_CLOSE
