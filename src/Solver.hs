{-# LANGUAGE TypeFamilies #-}

module Solver where

import Control.Monad.State
import Data.Functor.Rep
import Data.Kind (Type)
import Linear
import Linear.Covector

class Term term where
  type T term :: Type -> Type
  type S term :: Type -> Type
  type U term :: Type -> Type

  vf :: Num a => term a -> a -> T term a -> S term a
  control :: Num a => term a -> (a, a) -> U term a
  prod :: Num a => term a -> S term a -> U term a -> T term a

data SimpleODE v a = SimpleODE (a -> v a -> v a)

instance (Functor v) => Term (SimpleODE v) where
  type T (SimpleODE v) = v
  type S (SimpleODE v) = v
  type U (SimpleODE v) = V1

  vf (SimpleODE f) x t = f x t
  control _ (t0, t1) = V1 $ t1 - t0
  prod :: Num a => SimpleODE v a -> v a -> V1 a -> v a
  prod _ x (V1 dt) = dt *^ x

class Solver sol where
  init ::
    (Term term, Num a) =>
    term a ->
    (a, a) ->
    S term a ->
    sol (S term) a

  step ::
    (Term ode, Num a) =>
    ode a ->
    (T ode) a ->
    (a, a) ->
    StateT (sol (T ode) a) Maybe (T ode a)

data EulerSol v a = EulerSol

instance Solver EulerSol where
  init _ _ _ = EulerSol
  step ode y (t0, t1) = do
    let y' = vf ode t0 y
    let dt = control ode (t0, t1)
    return $ prod ode y' dt
