{-# LANGUAGE GADTs #-}
{-# LANGUAGE RecordWildCards #-}

module Solver.Class where

import Control.Monad.Identity (Identity)
import Control.Monad.RWS
import Control.Monad.State
import Data.Data (Proxy)
import Data.Functor.Rep
import Data.Kind (Type)
import Debug.Trace (traceShow, traceShowId)
import Linear
import Linear.Affine (Vector)
import Optics
import Term

class (Term de) => Solver sol de where
  initSolver ::
    (Num a) =>
    de a ->
    (a, a) ->
    T de a ->
    sol (T de) a

  step ::
    (Additive (T de), Floating a, Ord a, MonadState (sol (T de) a) m) =>
    de a ->
    (T de) a ->
    (a, a) ->
    m (T de a)

-- | Solvers which can offer an estimate of their error
class (Solver sol de) => ErrEst sol de where
  errorEstimate ::
    (Ord a, Floating a, Metric (T de)) =>
    a ->
    a ->
    de a ->
    sol (T de) a ->
    a
  errorOrder ::
    (Ord a, Floating a, Metric (T de)) =>
    de a ->
    sol (T de) a ->
    Int

class (Term ode, Solver solver ode) => Stepper stepper solver ode where
  passedDiscontinuity ::
    (Metric (T ode), Floating a) =>
    ode a ->
    stepper a ->
    solver (T ode) a ->
    Bool

  adaptStep ::
    (MonadState (stepper a) m, Ord a, (Metric (T ode), Floating a)) =>
    solver (T ode) a ->
    ode a ->
    (a, a) -> -- t0, t1
    (T ode a) -> -- y at t0
    (T ode a) -> -- y at t1 (candidate)

    -- | --- was candidate accepted?
    m (Maybe (T ode a), (a, a))

data ConstantStep a = ConstantStep {constantStepSize :: a}

instance (Solver solver de) => Stepper ConstantStep solver de where
  passedDiscontinuity _ _ _ = False
  adaptStep _ _ (_, t1) _y0 y1 = do
    dt <- gets constantStepSize
    return (Just y1, (t1, t1 + dt))

data SimpleAdaptStep a = SimpleAdaptStep
  { aTol :: a,
    rTol :: a
  }

instance
  (Solver solver de, ErrEst solver de) =>
  Stepper SimpleAdaptStep solver de
  where
  passedDiscontinuity _ _ _ = False

  adaptStep sol ode (t0, t1) _y0 y1 = do
    SimpleAdaptStep {..} <- get
    let q = errorOrder ode sol
    let err = (errorEstimate aTol rTol ode sol)
    let prop = max 0.25 . min 4 $ (t1 - t0) * (err) -- ^^ (-(q + 1))
    let h' = 0.9 * prop
    return $
      if err > 1
        then (Nothing, (t0, t0 + h'))
        else (Just y1, (t1, t1 + h'))

solveStep ::
  forall solver ode stepper a m.
  ( Solver solver ode,
    Stepper stepper solver ode,
    -- MonadState m,
    Floating a,
    Ord a,
    Show a
  ) =>
  ode a ->
  (T ode a) ->
  (a, a) ->
  State (stepper a, solver (T ode) a) ((a, a), Maybe (T ode a))
solveStep ode y0 (t0, t1) = do
  y1 <- zoom _2 $ step ode y0 (t0, t1)
  solver <- gets (view _2)
  (accepted, step') <- zoom _1 $ adaptStep solver ode (t0, t1) y0 y1
  return (step', accepted)

iterateWhileM :: (a -> Bool) -> (a -> State s a) -> a -> State s [a]
iterateWhileM pred f x =
  if pred x
    then (x :) <$> ((f x) >>= iterateWhileM pred f)
    else pure []

iterateSol :: (t, a) -> ((t, a) -> State s (t, Maybe a)) -> State s [(t, a)]
iterateSol (t, x) step = do
  (t', x') <- step (t, x)
  case x' of
    Nothing -> iterateSol (t', x) step
    Just x' -> ((t', x') :) <$> iterateSol (t', x') step

solve ::
  forall a v ode solver stepper.
  (Solver solver ode, Stepper stepper solver ode, Floating a, Ord a, Show a) =>
  ode a ->
  solver (T ode) a ->
  stepper a ->
  T ode a ->
  (a, a) ->
  [((a, a), T ode a)]
solve ode solver stepper y0 (ti, tf) =
  takeWhile (\((_, t), _) -> t <= tf) $
    evalState
      (iterateSol ((t0, t1), y0) go)
      (stepper, solver)
  where
    pred = (\((t, _), _) -> t <= tf)
    dt = (tf - ti) / 20
    (t0, t1) = (ti, ti + dt)
    go :: ((a, a), T ode a) -> State (stepper a, solver (T ode) a) ((a, a), Maybe (T ode a))
    go (step, y) = solveStep ode y step
