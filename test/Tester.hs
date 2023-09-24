{-# LANGUAGE ImpredicativeTypes #-}

module Tester where

import qualified Data.Vector as V
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
    forall sol stepper de a.
    (Solver sol de, Stepper stepper sol de) =>
    sol (T de) a ->
    stepper a ->
    SolParams sol de a ->
    Reference de a

data DETest where
  Compare ::
    forall de a.
    (Floating a, Show a, Ord a, Term de) =>
    [de a -> (a, a) -> T de a -> Interp de a] ->
    IVP de a ->
    -- | Tolerance
    a ->
    DETest
  DETest ::
    forall sol stepper ode a.
    ( Floating a,
      Show a,
      Ord a,
      Solver sol ode,
      Stepper stepper sol ode,
      Metric (T ode)
    ) =>
    -- | Solver to test
    sol (T ode) a ->
    -- | Parameters used by solver in test
    SolParams sol ode a ->
    stepper a ->
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
  run _ (DETest sol solParams st ivp tol) _ = do
    let IVP de exact (t0, t1) ts y0 = ivp
    let interps = solve de sol st solParams y0 (t0, t1)
    let exacts = exact <$> ts
    let ys = evalSol de ts interps
    let max_err = maximum . fmap norm $ zipWith (^-^) ys exacts
    return $
      if max_err / tol >= 1
        then testFailed ("unacceptable error " ++ show max_err)
        else testPassed ""
  testOptions = pure []
