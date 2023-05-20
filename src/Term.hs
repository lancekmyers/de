{-# LANGUAGE TypeFamilies #-}

module Term where

import Data.Kind (Type)
import Linear

class
  (Additive (T term), Additive (S term), Additive (U term)) =>
  Term term
  where
  type T term :: Type -> Type
  type S term :: Type -> Type
  type U term :: Type -> Type

  vf :: Floating a => term a -> a -> T term a -> S term a
  control :: Floating a => term a -> (a, a) -> U term a
  prod :: Floating a => term a -> S term a -> U term a -> T term a

data SimpleODE v a = SimpleODE (a -> v a -> v a)

instance (Additive v) => Term (SimpleODE v) where
  type T (SimpleODE v) = v
  type S (SimpleODE v) = v
  type U (SimpleODE v) = V1

  vf (SimpleODE f) x t = f x t
  control _ (t0, t1) = V1 $ t1 - t0
  prod :: Num a => SimpleODE v a -> v a -> V1 a -> v a
  prod _ x (V1 dt) = dt *^ x
