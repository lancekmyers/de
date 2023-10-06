{-# LANGUAGE GADTs #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}

module Solver.Class
  ( StepperPID,
    basicI,
    pi42,
    pi33,
    pi34,
    h211PI,
    h312PID,
    ConstantStepper (..),
    ErrEst (..),
    Solver (..),
    Stepper (..),
    solve,
  )
where

import Control.Monad.Identity (Identity)
import Control.Monad.RWS
import Control.Monad.Reader (ReaderT)
import Control.Monad.State
import Data.Data (Proxy)
import Data.Functor.Rep
import Data.Kind (Type)
import Debug.Trace (traceShow, traceShowId)
import GHC.Generics (Generic)
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
    (MonadState (stepper a) m, Ord a, (Metric (T ode), Floating a, Real a)) =>
    SolState solver ode a ->
    solver (T ode) a ->
    stepper a ->
    ode a ->
    (a, a) -> -- t0, t1
    --- (T ode a) -> -- y at t0
    --- (T ode a) -> -- y at t1 (candidate)

    -- | --- was candidate accepted?
    m (Either (a, a) (a, a))

data ConstantStepper a = ConstantStepper

instance
  (Solver solver de) =>
  Stepper ConstantStepper solver de
  where
  adaptStep solState sol _ ode (t0, t1) = do
    let h = t1 - t0
    return $ Right (t1, t1 + h)

data StepperPID a = StepperPID
  { prevErr :: a,
    prevPrevErr :: a,
    beta_1 :: a,
    beta_2 :: a,
    beta_3 :: a
  }
  deriving (Generic)

mkStepperPID :: Floating a => a -> a -> a -> StepperPID a
mkStepperPID p i d =
  StepperPID
    { prevErr = 1,
      prevPrevErr = 1,
      beta_1 = beta_1,
      beta_2 = beta_3,
      beta_3 = beta_3
    }
  where
    beta_1 = p + i + d
    beta_2 = -(p + 2 * d)
    beta_3 = d

mkStepperPI :: Floating a => a -> a -> StepperPID a
mkStepperPI p i = mkStepperPID p i 0

mkStepperI :: Floating a => a -> StepperPID a
mkStepperI i = mkStepperPID 0 i 0

defaultStepperPI :: Floating a => StepperPID a
defaultStepperPI = mkStepperPI 0.2 1.0

instance
  (Solver solver de, ErrEst solver de) =>
  Stepper StepperPID solver de
  where
  adaptStep solState sol _ ode (t0, t1) = do
    StepperPID {..} <- get
    let q = errorOrder ode sol
    let err = (errorEstimate ode sol solState)
    let prop = err ** beta_1 * prevErr ** beta_2 * prevPrevErr ** beta_3
    -- limit step factor
    let prop' = 1 + atan (prop - 1)
    let h' = prop' * (t1 - t0)

    if err > 1
      then return $ Left (t0, t0 + h')
      else do
        assign #prevPrevErr prevErr
        assign #prevErr err
        return $ Right (t1, t1 + h')

basicI :: Floating a => StepperPID a
basicI =
  StepperPID
    { prevErr = 1,
      prevPrevErr = 1,
      beta_1 = 1.0,
      beta_2 = 0.0,
      beta_3 = 0.0
    }

pi42 :: Floating a => StepperPID a
pi42 =
  StepperPID
    { prevErr = 1,
      prevPrevErr = 1,
      beta_1 = 0.6,
      beta_2 = -0.2,
      beta_3 = 0.0
    }

pi33 :: Floating a => StepperPID a
pi33 =
  StepperPID
    { prevErr = 1,
      prevPrevErr = 1,
      beta_1 = 2 / 3,
      beta_2 = -1 / 3,
      beta_3 = 0.0
    }

pi34 :: Floating a => StepperPID a
pi34 =
  StepperPID
    { prevErr = 1,
      prevPrevErr = 1,
      beta_1 = 0.7,
      beta_2 = -0.4,
      beta_3 = 0.0
    }

h211PI :: Floating a => StepperPID a
h211PI =
  StepperPID
    { prevErr = 1,
      prevPrevErr = 1,
      beta_1 = 1 / 6,
      beta_2 = 1 / 6,
      beta_3 = 0.0
    }

h312PID :: Floating a => StepperPID a
h312PID =
  StepperPID
    { prevErr = 1,
      prevPrevErr = 1,
      beta_1 = 1 / 18,
      beta_2 = 1 / 9,
      beta_3 = 1 / 18
    }

solveStep ::
  forall solver ode stepper a.
  ( Solver solver ode,
    Stepper stepper solver ode,
    Floating a,
    Real a,
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
  (Solver solver ode, Stepper stepper solver ode, Real a, Floating a, Ord a, Show a) =>
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
