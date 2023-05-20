{-# LANGUAGE DataKinds #-}
{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module Solver where

import Control.Applicative
import Control.Monad
import Control.Monad.Identity (Identity)
import Control.Monad.State
import Data.Functor.Rep
import Data.Kind (Type)
import GHC.TypeLits (Nat, type (+))
import Linear
import Linear.Covector
import Linear.V
import Linear.Vector
import Optics

class Term term where
  type T term :: Type -> Type
  type S term :: Type -> Type
  type U term :: Type -> Type

  vf :: Num a => term a -> a -> T term a -> S term a
  control :: Num a => term a -> (a, a) -> U term a
  prod :: Num a => term a -> S term a -> U term a -> T term a

data SimpleODE v a = SimpleODE (a -> v a -> v a)

instance (Additive v) => Term (SimpleODE v) where
  type T (SimpleODE v) = v
  type S (SimpleODE v) = v
  type U (SimpleODE v) = V1

  vf (SimpleODE f) x t = f x t
  control _ (t0, t1) = V1 $ t1 - t0
  prod :: Num a => SimpleODE v a -> v a -> V1 a -> v a
  prod _ x (V1 dt) = dt *^ x

class (Term de) => Solver sol de where
  initSolver ::
    (Num a) =>
    de a ->
    (a, a) ->
    T de a ->
    sol (T de) a

  step ::
    (Additive (T de), Num a, MonadState (sol (T de) a) m) =>
    de a ->
    (T de) a ->
    (a, a) ->
    m (T de a)

class (Solver sol de) => ErrEst sol de where
  errorEstimate :: de a -> sol (T de) a -> T de a

data EulerSol v a = EulerSol

instance (Term ode) => Solver EulerSol ode where
  initSolver _ _ _ = EulerSol
  step ode y (t0, t1) = do
    let y' = vf ode t0 y
    let dt = control ode (t0, t1)
    return $ y ^+^ prod ode y' dt

data Heun v a = Heun

class (Term ode, Solver solver ode) => Stepper stepper solver ode where
  passedDiscontinuity :: ode a -> stepper a -> solver (T ode) a -> Bool

  initStepper ::
    solver (T ode) a ->
    ode a ->
    (a, a) ->
    (T ode a) ->
    a ->
    stepper a

  adaptStep ::
    (MonadState (stepper a) m, Num a) =>
    solver (T ode) a ->
    ode a ->
    (a, a) -> -- t0, t1
    (T ode a) -> -- y at t0
    (T ode a) -> -- y at t1 (candidate)

    -- | --- was candidate accepted?
    m (Bool, (a, a))

data ConstantStep a = ConstantStep {constantStepSize :: a}

instance (Solver solver de) => Stepper ConstantStep solver de where
  passedDiscontinuity _ _ _ = False
  initStepper _ _ _ _ dt = ConstantStep dt
  adaptStep _ _ (_, t1) _ _ = do
    dt <- gets constantStepSize
    return (True, (t1, t1 + dt))

solveStep ::
  forall solver ode stepper a m.
  ( Solver solver ode,
    Additive (T ode),
    Stepper stepper solver ode,
    -- MonadState m,
    Num a
  ) =>
  ode a ->
  (T ode a) ->
  (a, a) ->
  State (stepper a, solver (T ode) a) ((a, a), T ode a)
solveStep ode y0 (t0, t1) = do
  y1 <- zoom _2 $ step ode y0 (t0, t1)
  solver <- gets (view _2)
  (accepted, step') <- zoom _1 $ adaptStep solver ode (t0, t1) y0 y1
  return (step', y1)

iterateWhileM :: (a -> Bool) -> (a -> State s a) -> a -> State s [a]
iterateWhileM pred f x =
  if pred x
    then (x :) <$> ((f x) >>= iterateWhileM pred f)
    else pure []

solve ::
  forall a v.
  (Additive v, Floating a, Ord a) =>
  SimpleODE v a ->
  v a ->
  (a, a) ->
  [((a, a), v a)]
solve ode y0 (ti, tf) = fst $ runState (iterateWhileM pred go ((t0, t1), y0)) (stepper, solver)
  where
    pred = (\((t, _), _) -> t <= tf)
    dt = (tf - ti) / 20
    (t0, t1) = (ti, ti + dt)
    stepper = initStepper solver ode (t0, t1) y0 dt
    solver = initSolver ode (t0, t1) y0
    go :: ((a, a), v a) -> State (ConstantStep a, EulerSol v a) ((a, a), v a)
    go (step, y) = solveStep ode y step
