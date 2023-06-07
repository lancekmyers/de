module Tester where

import Interpolate
import Linear
import Solver
import Term
import Test.Tasty.Options (OptionDescription)
import Test.Tasty.Providers

data IVP de a = IVP
  { de :: de a,
    exact :: a -> T de a,
    span :: (a, a),
    interiorPts :: [a],
    y0 :: T de a
  }

data DETest where
  DETest ::
    forall sol stepper ode a.
    ( Floating a,
      Show a,
      Ord a,
      Solver sol ode,
      Stepper stepper sol ode,
      Metric (T ode)
    ) =>
    sol (T ode) a ->
    SolParams sol ode a ->
    stepper a ->
    IVP ode a ->
    a ->
    DETest

instance IsTest DETest where
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
