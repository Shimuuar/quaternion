{-# LANGUAGE ViewPatterns          #-}
{-# LANGUAGE FlexibleInstances     #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE TypeFamilies          #-}
-- |
-- Quaternions
module Data.Quaternion (
    QuaternionG(..)
  , Quaternion
  , scaleQuat
  , conjugate
  , realQuat
  , quaternionNorm
  , quaternionNormSq
  ) where

import qualified Data.Vector.Fixed         as F
import           Data.Vector.Fixed         (Dim,Vector,N4)
import           Data.Vector.Fixed.Unboxed (Vec4)



-- | Wrapper new for any 4-element vector.
newtype QuaternionG v a = QuaternionG (v a)

-- | Quaternions wrapping unboxed vectors.
type Quaternion = QuaternionG Vec4


type instance Dim (QuaternionG v) = N4

instance (Vector v a, Dim v ~ N4) => Vector (QuaternionG v) a where
  inspect (QuaternionG v) f = F.inspect v f
  construct = fmap QuaternionG F.construct
  {-# INLINE inspect #-}
  {-# INLINE construct #-}

instance (Vector v a, Dim v ~ N4, Floating a) => Num (QuaternionG v a) where
  (+) = F.zipWith (+)
  (-) = F.zipWith (-)
  (F.convert -> (w1,x1,y1,z1)) * (F.convert -> (w2,x2,y2,z2))
    = F.mk4 (w1*w2 - x1*x2 - y1*y2 - z1*z2)
            (w1*x2 + x1*w2 + y1*z2 - z1*y2)
            (w1*y2 - x1*z2 + y1*w2 + z1*x2)
            (w1*z2 + x1*y2 - y1*x2 + z1*w2)
  {-# INLINE (+) #-}
  {-# INLINE (-) #-}
  {-# INLINE (*) #-}

  abs    = realQuat . quaternionNorm
  signum = realQuat . signum . F.head
  {-# INLINE abs #-}
  {-# INLINE signum #-}

  fromInteger i = realQuat (fromInteger i)
  {-# INLINE fromInteger #-}

instance (Vector v a, Dim v ~ N4, Floating a) => Fractional (QuaternionG v a) where
  v / u   = recip (quaternionNormSq v) `scaleQuat` v * conjugate u
  recip v = recip (quaternionNormSq v) `scaleQuat` conjugate v
  {-# INLINE (/)   #-}
  {-# INLINE recip #-}
  fromRational = realQuat . fromRational
  {-# INLINE fromRational #-}


-- | Conjugate quaternion
conjugate :: (Vector v a, Dim v ~ N4, Num a) => QuaternionG v a -> QuaternionG v a
conjugate = F.imap (\i a -> if i == 0 then a else negate a)
{-# INLINE conjugate #-}

-- | Scale quaternion
scaleQuat :: (Vector v a, Dim v ~ N4, Num a) => a -> QuaternionG v a -> QuaternionG v a
scaleQuat a = F.map (* a)
{-# INLINE scaleQuat #-}

-- | Construct quaternion from real part
realQuat :: (Vector v a, Dim v ~ N4, Num a) => a -> QuaternionG v a
realQuat x = F.mk4 x 0 0 0
{-# INLINE realQuat #-}

-- | Norm of quaternion
quaternionNorm :: (Vector v a, Dim v ~ N4, Floating a) => QuaternionG v a -> a
quaternionNorm = sqrt . quaternionNormSq
{-# INLINE quaternionNorm #-}

-- | Square of quaternion's norm
quaternionNormSq :: (Vector v a, Dim v ~ N4, Num a) => QuaternionG v a -> a
quaternionNormSq = F.sum . F.map (\x -> x*x)
{-# INLINE quaternionNormSq #-}
