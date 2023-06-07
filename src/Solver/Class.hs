{-# LANGUAGE GADTs #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}

module Solver.Class where

import Control.Monad.Identity (Identity)
import Control.Monad.RWS
import Control.Monad.Reader (ReaderT)
import Control.Monad.State
import Data.Data (Proxy)
import Data.Functor.Rep
import Data.Kind (Type)
import Debug.Trace (traceShow, traceShowId)
import Interpolate
import Linear
import Linear.Affine (Vector)
import Optics
import Term

class (Term de) => Solver sol de where
  type SolState sol de :: Type -> Type
  type SolParams sol de :: Type -> Type

  initSolver ::
    (Num a) =>
    de a ->
    (a, a) ->
    T de a ->
    sol (T de) a ->
    SolParams sol de a ->
    SolState sol de a

  step ::
    ( Additive (T de),
      Floating a,
      Ord a,
      MonadState (SolState sol de a) m,
      MonadReader (SolParams sol de a) m
    ) =>
    sol (T de) a ->
    de a ->
    (T de) a ->
    (a, a) ->
    m (Interp de a)

-- | Solvers which can offer an estimate of their error
class (Solver sol de) => ErrEst sol de where
  errorEstimate ::
    (Ord a, Floating a, Metric (T de)) =>
    de a ->
    sol (T de) a ->
    SolState sol de a ->
    a
  errorOrder ::
    (Ord a, Floating a, Metric (T de)) =>
    de a ->
    sol (T de) a ->
    Int

class (Term ode, Solver solver ode) => Stepper stepper solver ode where
  adaptStep ::
    (MonadState (stepper a) m, Ord a, (Metric (T ode), Floating a)) =>
    SolState solver ode a ->
    solver (T ode) a ->
    stepper a ->
    ode a ->
    (a, a) -> -- t0, t1
    --- (T ode a) -> -- y at t0
    --- (T ode a) -> -- y at t1 (candidate)

    -- | --- was candidate accepted?
    m (Either (a, a) (a, a))

data ConstantStep a = ConstantStep {constantStepSize :: a}

instance (Solver solver de) => Stepper ConstantStep solver de where
  -- passedDiscontinuity _ _ _ = False
  adaptStep _ _ _ __ (t0, t1) = do
    dt <- gets constantStepSize
    return (Right (t1, t1 + dt))

-- | I controller
data SimpleAdaptStep a = SimpleAdaptStep

instance
  (Solver solver de, ErrEst solver de) =>
  Stepper SimpleAdaptStep solver de
  where
  -- passedDiscontinuity _ _ _ = False

  adaptStep solState sol _ ode (t0, t1) = do
    let q = errorOrder ode sol
    let err = (errorEstimate ode sol solState)
    let prop = max 0.25 . min 4 $ (t1 - t0) * (err) ^^ (-(q + 1))
    let h' = 0.9 * prop * (t1 - t0)
    return $
      if err > 1
        then Left (t0, t0 + h')
        else Right (t1, t1 + h')

solveStep ::
  forall solver ode stepper a.
  ( Solver solver ode,
    Stepper stepper solver ode,
    Floating a,
    Ord a,
    Show a
  ) =>
  ode a ->
  solver (T ode) a ->
  stepper a ->
  (T ode a) ->
  (a, a) ->
  RWS
    (SolParams solver ode a)
    [Interp ode a]
    (stepper a, SolState solver ode a)
    ( ((a, a), T ode a)
    )
solveStep ode sol stepper y0 (t0, t1) = do
  -- now returns Interp, not y
  interpolant <- zoom _2 $ step sol ode y0 (t0, t1)

  solState <- gets (view _2)
  step <-
    zoom _1 $
      adaptStep solState sol stepper ode (t0, t1)

  case step of
    Left step' -> solveStep ode sol stepper y0 step'
    Right step' -> tell [interpolant] >> return (step', rightMost interpolant)

iterateWhileM :: Monad m => (a -> Bool) -> (a -> m a) -> a -> m a
iterateWhileM pred f x = do
  x' <- f x
  if pred x'
    then iterateWhileM pred f x'
    else return x'

-- iterateSol :: (t, a) -> ((t, a) -> State s (t, Maybe a)) -> State s [(t, a)]
iterateSol :: Monad m => (a, b) -> ((a, b) -> m (a, Maybe b)) -> m [(a, b)]
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
  SolParams solver ode a ->
  T ode a ->
  (a, a) ->
  [Interp ode a] -- [((a, a), T ode a)]
solve ode solver stepper solParams y0 (ti, tf) =
  snd $
    execRWS
      (iterateWhileM pred go ((t0, t1), y0))
      solParams
      (stepper, solState)
  where
    pred = (\((t, _), _) -> t <= tf)
    dt = (tf - ti) / 20
    (t0, t1) = (ti, ti + dt)

    solState = initSolver ode (ti, tf) y0 solver solParams
    go ::
      ((a, a), T ode a) ->
      RWS
        (SolParams solver ode a)
        [Interp ode a]
        (stepper a, SolState solver ode a)
        ((a, a), T ode a)
    go (step, y) = solveStep ode solver stepper y step
