{-# LANGUAGE ImpredicativeTypes #-}

module Tester where

import Data.Vector qualified as V
import Interpolate
import Linear
import Solver
import Term
import Test.Tasty.Options (OptionDescription)
import Test.Tasty.Providers

linspace :: (Ord a, Floating a) => Int -> a -> a -> [a]
linspace n a b
  | n < 1 = []
  | b < a = []
  | otherwise = take n $ iterate (+ w) a
  where
    w = (b - a) / fromIntegral n

data IVP de a = IVP
  { de :: de a,
    exact :: a -> T de a,
    span :: (a, a),
    interiorPts :: [a],
    y0 :: T de a
  }

data Reference de a where
  Exact :: (a -> T de a) -> Reference de a
  Numerical ::
    (Term de) =>
    StepIntegrator (T de) a ->
    StepController a ->
    Reference de a

data DETest where
  Compare ::
    forall de a.
    (Floating a, Show a, Ord a, Term de) =>
    [de a -> (a, a) -> T de a -> Interp (T de) a] ->
    IVP de a ->
    -- | Tolerance
    a ->
    DETest
  DETest ::
    forall sol stepper ode a.
    ( Floating a,
      Real a,
      Show a,
      Ord a,
      Term ode,
      Metric (T ode)
    ) =>
    -- | Solver to test
    (ode a -> StepIntegrator (T ode) a) ->
    StepController a ->
    IVP ode a ->
    -- | Tolerance
    a ->
    DETest

{-
  Also need a way to generate html + plots for these
-}

instance IsTest DETest where
  run _ (Compare sols ivp tol) _ = do
    let IVP de exact (t0, t1) ts y0 = ivp
    let solutions = [solver de (t0, t1) y0 | solver <- sols]
    let pts = V.fromList $ linspace 10 t0 t1
    let solvedPts = [V.map (\t -> interp de t solution) pts | solution <- solutions]
    let go (xs, ys) = V.maximum $ V.map norm $ V.zipWith (^-^) xs ys
    let max_err = maximum $ go <$> ((,) <$> solvedPts <*> solvedPts)
    return $
      if max_err >= tol
        then testFailed ("unacceptable error " ++ show max_err)
        else testPassed ""
  run _ (DETest stepInt timeCont ivp tol) _ = do
    let IVP de exact (t0, t1) ts y0 = ivp
    let interps = runIntegration (stepInt de) timeCont (y0, TimeStep {t = t0, delta = (t1 - t0)}) (last ts)
    let exacts = exact <$> ts
    let ys = evalSol de ts interps
    let max_err = maximum . fmap norm $ zipWith (^-^) ys exacts
    let last_err = norm $ last ys ^-^ last exacts
    return $
      if max_err / tol >= 1
        then testFailed ("unacceptable error " ++ show max_err ++ " | " ++ show last_err)
        else testPassed ""
  testOptions = pure []
