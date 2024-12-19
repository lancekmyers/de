{-# LANGUAGE GADTs #-}
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
    type StepController,
    type StepIntegrator,
    ErrorEstimate (..),
    TimeStep (..),
    runIntegration,
    solvingMachine,
  )
where

import Control.Arrow (Arrow (..))
import Control.Category (Category (..), (>>>))
import Control.Monad.Identity (Identity)
import Control.Monad.RWS
import Control.Monad.Reader (ReaderT)
import Control.Monad.State
import Data.Data (Proxy)
import Data.Functor.Rep
import Data.Kind (Type)
import Data.Machine (Is, Machine, Plan, PlanT, auto, await, construct, run, source, taking, yield)
import Data.Machine.Mealy (Mealy (..), unfoldMealy)
import Data.Machine.Process (Process, takingWhile, (~>))
import Debug.Trace (traceShow, traceShowId)
import GHC.Generics (Generic)
import Interpolate
import Linear
import Optics hiding (Is)
import Term
import Prelude hiding ((.))

runIntegration ::
  forall a v solver stepper.
  (Real a, Floating a, Ord a, Show a, Additive v) =>
  StepIntegrator v a ->
  StepController a ->
  (v a, TimeStep a) ->
  a ->
  [Interp v a]
runIntegration stepIntegrator stepController (y0, t0) tf =
  run $
    source [(y0, t0)] ~> solveMach ~> takingWhile (\(Poly (_t0, t1) _) -> tf > t1)
  where
    solveMach = solvingMachine stepIntegrator stepController

buildStep ::
  StepIntegrator v a ->
  StepController a ->
  Mealy (v a, TimeStep a) (Either (TimeStep a) (Interp v a, TimeStep a))
buildStep sol stp = fmap sequence $ (sol &&& arr snd) >>> lAssoc >>> second stp
  where
    lAssoc = arr $ \((int, err), tst) -> (int, (err, tst))

solvingMachine ::
  forall a v k.
  (Show a, Num a, Additive v) =>
  StepIntegrator v a ->
  StepController a ->
  Machine (Is (v a, TimeStep a)) (Interp v a) -- Plan (Is (v a, TimeStep a)) (Interp v a) ()
solvingMachine sol stp = construct $ await >>= loop mealy
  where
    lAssoc = arr $ \((int, err), tst) -> (int, (traceShow err err, tst))
    -- integrate step, then check error and reject/accept with new step
    mealy = (sol &&& arr snd) >>> lAssoc >>> second stp >>> arr sequence
    loop ::
      forall k m.
      (Monad m) =>
      Mealy (v a, TimeStep a) (Either (TimeStep a) (Interp v a, TimeStep a)) ->
      (v a, TimeStep a) ->
      PlanT (Is (v a, TimeStep a)) (Interp v a) m ()
    loop mealy (y0, h) = do
      let (ret, mealy') = runMealy mealy (y0, h)
      case traceShow h ret of
        -- accepted
        Right (interp, h') -> yield interp >> loop mealy' (rightMost interp, traceShow "accept" h')
        -- rejected
        Left h' -> loop mealy' (y0, traceShow "reject" h')

interpTimeStep :: (Num a) => Interp v a -> TimeStep a
interpTimeStep (Poly (t0, t1) _) = TimeStep {t = t0, delta = t1 - t0}

type StepIntegrator v a = Mealy (v a, TimeStep a) (Interp v a, ErrorEstimate a)

data ErrorEstimate a = ErrorEstimate a Int
  deriving (Show)

data TimeStep a = TimeStep {t :: a, delta :: a}
  deriving (Show)

type StepController a = Mealy (ErrorEstimate a, TimeStep a) (Either (TimeStep a) (TimeStep a))

constantStepper :: (Num a) => StepController a
constantStepper = unfoldMealy go ()
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

mkStepperPID :: (Floating a, Ord a) => a -> a -> a -> StepController a
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

mkStepperPI :: (Floating a, Ord a) => a -> a -> StepController a
mkStepperPI p i = mkStepperPID p i 0

mkStepperI :: (Floating a, Ord a) => a -> StepController a
mkStepperI i = mkStepperPID 0 i 0

defaultStepperPI :: (Floating a, Ord a) => StepController a
defaultStepperPI = mkStepperPI 0.2 1.0

pidStepController :: forall a. (Ord a, Floating a) => StepperPID a -> StepController a
pidStepController pid = unfoldMealy go pid
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

basicI :: (Floating a, Ord a) => StepController a
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

pi42 :: (Floating a, Ord a) => StepController a
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

pi33 :: (Floating a, Ord a) => StepController a
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

pi34 :: (Floating a, Ord a) => StepController a
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

h211PI :: (Floating a, Ord a) => StepController a
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

h312PID :: (Floating a, Ord a) => StepController a
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
