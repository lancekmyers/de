{-# LANGUAGE TypeFamilies #-}

module Term (Term (..), SimpleODE (..), JacODE (..)) where

import Data.Kind (Type)
import Linear

class
  (Metric (T term), Metric (S term), Metric (U term)) =>
  Term term
  where
  type T term :: Type -> Type
  type S term :: Type -> Type
  type U term :: Type -> Type

  vf :: (Floating a) => term a -> a -> T term a -> S term a
  control :: (Floating a) => term a -> (a, a) -> U term a
  prod :: (Floating a) => term a -> S term a -> U term a -> T term a

newtype SimpleODE v a = SimpleODE (a -> v a -> v a)

instance (Metric v) => Term (SimpleODE v) where
  type T (SimpleODE v) = v
  type S (SimpleODE v) = v
  type U (SimpleODE v) = V1
  control _ (t0, t1) = V1 $ t1 - t0
  prod :: (Num a) => SimpleODE v a -> v a -> V1 a -> v a
  prod _ x (V1 dt) = dt *^ x

  vf (SimpleODE f) = f

-- \| ODE that comes with jacobian
data JacODE v a = JacODE
  { jac :: a -> v a -> v (v a),
    runvf :: a -> v a -> v a
  }

instance (Additive v, Metric v) => Term (JacODE v) where
  type S (JacODE v) = v
  type T (JacODE v) = v
  type U (JacODE v) = V1
  vf ::
    (Additive v, Floating a) =>
    JacODE v a ->
    a ->
    T (JacODE v) a ->
    S (JacODE v) a
  vf JacODE {jac, runvf} = runvf
  control ::
    (Additive v, Floating a) =>
    JacODE v a ->
    (a, a) ->
    U (JacODE v) a
  control _ (t1, t2) = V1 (t2 - t1)

  prod ::
    (Additive v, Floating a) =>
    JacODE v a ->
    S (JacODE v) a ->
    U (JacODE v) a ->
    T (JacODE v) a
  prod _ x (V1 dt) = dt *^ x
