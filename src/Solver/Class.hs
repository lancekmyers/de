{-# LANGUAGE DerivingVia #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE LambdaCase #-}
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
    constantStepper,
    controlledStep,
    type StepController,
    type Solver,
    SolverInfo (..),
    SolverErr (..),
    ErrorEstimate (..),
    TimeStep (..),
    runIntegration,
  )
where

import Control.Arrow (Arrow (..))
import Control.Category (Category (..), (>>>))
import Control.Monad.Except (Except, ExceptT)
import Control.Monad.Identity (Identity)
import Control.Monad.RWS
import Control.Monad.Reader (ReaderT)
import Control.Monad.State
import Control.Monad.Writer.Lazy (WriterT)
import Data.Data (Proxy)
import Data.Functor.Compose (Compose (..))
import Data.Functor.Rep
import Data.Kind (Type)
import Data.Machine (Is, Machine, MachineT (runMachineT), Plan, PlanT, ProcessT, Step (..), auto, await, construct, run, source, supply, taking, yield)
import Data.Machine.Mealy (Mealy (..), unfoldMealy)
import Data.Machine.MealyT (MealyT (runMealyT), upgrade)
import Data.Machine.Process (Process, takingWhile, (~>))
import Data.Machine.Runner (runT)
import Data.Semigroup (Sum)
import Debug.Trace (traceShow, traceShowId)
import GHC.Generics (Generic, Generically)
import Interpolate
import Linear
import Optics hiding (Is)
import Term
import Prelude hiding ((.))

runIntegration ::
  forall a v m.
  (Real a, Floating a, Ord a, Show a, Additive v, Monad m) =>
  Solver () m v a ->
  (v a, a) ->
  a ->
  a ->
  WriterT SolverInfo (ExceptT SolverErr m) [Interp v a]
runIntegration solver (y0, t0) h tf =
  runT $
    source [(y0, TimeStep t0 h)] ~> (snd <$> solver) ~> takingWhile (\(Poly (_t0, t1) _) -> tf >= t1)

controlledStep ::
  forall a v.
  (Show a, Num a, Additive v) =>
  StepController Identity a ->
  Solver (ErrorEstimate a) Identity v a ->
  Solver () Identity v a
controlledStep stp sol = construct (await >>= go sol stp)
  where
    go ::
      Solver (ErrorEstimate a) Identity v a ->
      StepController Identity a ->
      (v a, TimeStep a) ->
      PlanT (Is (v a, TimeStep a)) ((), Interp v a) (WriterT SolverInfo (Except SolverErr)) ()
    go sol stpMeal (y0, t) = do
      ((err, interp), sol') <- lift (runMachineT (supply [(y0, t)] sol)) >>= (\(Yield x k) -> pure (x, k))

      (t', stpMeal') <- lift $ runMealyT stp (err, t)

      case t' of
        Left t' -> go sol' stpMeal' (y0, t')
        Right t' -> yield ((), interp)

interpTimeStep :: (Num a) => Interp v a -> TimeStep a
interpTimeStep (Poly (t0, t1) _) = TimeStep {t = t0, delta = t1 - t0}

-- | A solver takes in a point and a time step and returns an interpolant over the interval,
-- along with some auxiliary information.
-- The solver is allowed to perform effects, in particular log information
-- about the solving process and throw errors.
type Solver i m v a =
  ProcessT (WriterT SolverInfo (ExceptT SolverErr m)) (v a, TimeStep a) (i, Interp v a)

-- | Information about a particular step taken by an integrator.
-- No information, but may later contain
--   - rejection count
--   - truncation error estimates
--   - convergence rate for newton solvers
data StepInfo = StepInfo
  deriving (Show)

data SolverInfo = SolverInfo Int [StepInfo]
  deriving (Show)

instance Semigroup SolverInfo where
  (<>) :: SolverInfo -> SolverInfo -> SolverInfo
  (SolverInfo n1 s1) <> (SolverInfo n2 s2) =
    SolverInfo (n1 + n2) (s1 <> s2)

instance Monoid SolverInfo where
  mempty :: SolverInfo
  mempty = SolverInfo 0 []

-- | Error thrown by solver when an issue occurs.
-- For now just a plain string, but will later include more info.
-- Could be that the newton solver fails to converge or that the stepper
-- has too many rejections.
data SolverErr = SolverErr String

data ErrorEstimate a = ErrorEstimate a Int
  deriving (Show)

data TimeStep a = TimeStep {t :: a, delta :: a}
  deriving (Show)

-- | Step size controller
-- This will adjust step size based on estimated error.
-- The returned value indicates acceptance or rejection of the proposed
-- timestep and a new time step to try.
type StepController m a =
  MealyT
    (WriterT SolverInfo (ExceptT SolverErr m))
    (ErrorEstimate a, TimeStep a)
    (Either (TimeStep a) (TimeStep a))

constantStepper :: (Num a) => StepController Identity a
constantStepper = upgrade $ unfoldMealy go ()
  where
    go _ (_errEst, TimeStep {t, delta}) =
      let t' = t + delta
          step = TimeStep {t = t', delta}
       in (Right step, ())

data StepperPID a = StepperPID
  { prevErr :: a,
    prevPrevErr :: a,
    rejections :: Int,
    beta_1 :: a,
    beta_2 :: a,
    beta_3 :: a
  }
  deriving (Generic)

mkStepperPID :: (Floating a, Ord a) => a -> a -> a -> StepController Identity a
mkStepperPID p i d =
  pidStepController
    StepperPID
      { prevErr = 1,
        prevPrevErr = 1,
        rejections = 0,
        beta_1 = beta_1,
        beta_2 = beta_3,
        beta_3 = beta_3
      }
  where
    beta_1 = p + i + d
    beta_2 = -(p + 2 * d)
    beta_3 = d

mkStepperPI :: (Floating a, Ord a) => a -> a -> StepController Identity a
mkStepperPI p i = mkStepperPID p i 0

mkStepperI :: (Floating a, Ord a) => a -> StepController Identity a
mkStepperI i = mkStepperPID 0 i 0

defaultStepperPI :: (Floating a, Ord a) => StepController Identity a
defaultStepperPI = mkStepperPI 0.2 1.0

pidStepController :: forall a. (Ord a, Floating a) => StepperPID a -> StepController Identity a
pidStepController pid = upgrade $ unfoldMealy go pid
  where
    go :: StepperPID a -> (ErrorEstimate a, TimeStep a) -> (Either (TimeStep a) (TimeStep a), StepperPID a)
    go pid@StepperPID {..} (ErrorEstimate err _, TimeStep {..})
      | rejections > 5 = error "too many rejections"
      -- \| err < 1e-4 = error "error is too small"
      | err > 1 = (Left (TimeStep {t = t, delta = delta'}), pid {rejections = rejections + 1})
      | otherwise =
          ( Right (TimeStep {t = t + delta, delta = delta'}),
            pid {rejections = 0, prevPrevErr = prevErr, prevErr = 1e-2 + err}
          )
      where
        prop = ((err) ** beta_1) * (prevErr ** beta_2)
        -- \* (prevPrevErr ** beta_3)
        -- limit step factor
        prop' = min 2 . max 0.5 $ 1 / prop
        delta' = prop' * delta

basicI :: (Floating a, Ord a) => StepController Identity a
basicI =
  pidStepController
    StepperPID
      { prevErr = 1,
        prevPrevErr = 1,
        rejections = 0,
        beta_1 = 1.0,
        beta_2 = 0.0,
        beta_3 = 0.0
      }

pi42 :: (Floating a, Ord a) => StepController Identity a
pi42 =
  pidStepController
    StepperPID
      { prevErr = 1,
        prevPrevErr = 1,
        rejections = 0,
        beta_1 = 0.6,
        beta_2 = -0.2,
        beta_3 = 0.0
      }

pi33 :: (Floating a, Ord a) => StepController Identity a
pi33 =
  pidStepController
    StepperPID
      { prevErr = 1,
        prevPrevErr = 1,
        rejections = 0,
        beta_1 = 2 / 3,
        beta_2 = -1 / 3,
        beta_3 = 0.0
      }

pi34 :: (Floating a, Ord a) => StepController Identity a
pi34 =
  pidStepController
    StepperPID
      { prevErr = 1,
        prevPrevErr = 1,
        rejections = 0,
        beta_1 = 0.7,
        beta_2 = -0.4,
        beta_3 = 0.0
      }

h211PI :: (Floating a, Ord a) => StepController Identity a
h211PI =
  pidStepController
    StepperPID
      { prevErr = 1,
        prevPrevErr = 1,
        rejections = 0,
        beta_1 = 1 / 6,
        beta_2 = 1 / 6,
        beta_3 = 0.0
      }

h312PID :: (Floating a, Ord a) => StepController Identity a
h312PID =
  pidStepController
    StepperPID
      { prevErr = 1,
        prevPrevErr = 1,
        rejections = 0,
        beta_1 = 1 / 18,
        beta_2 = 1 / 9,
        beta_3 = 1 / 18
      }
