{-# LANGUAGE ImpredicativeTypes #-}

module Tester where

import Control.Monad.Except (runExcept)
import Control.Monad.Identity (Identity)
import Control.Monad.Writer.Lazy (runWriterT)
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
    Solver () Identity (T de) a ->
    Reference de a

data DETest where
  DETest ::
    forall sol stepper ode a i.
    ( Floating a,
      Real a,
      Show a,
      Ord a,
      Term ode,
      Metric (T ode)
    ) =>
    -- | Solver to test
    (ode a -> Solver i Identity (T ode) a) ->
    IVP ode a ->
    -- | Initial time step
    a ->
    -- | Tolerance
    a ->
    DETest

{-
  Also need a way to generate html + plots for these
-}

instance IsTest DETest where
  run _ (DETest sol ivp h tol) _ = do
    let IVP de exact (t0, tf) ts y0 = ivp
    let interps' = runIntegration ((\(_, a) -> ((), a)) <$> sol de) (y0, t0) h tf
    interps <- case runExcept (runWriterT interps') of
      Left (SolverErr err) -> error err
      Right (x, _info) -> pure x
    let exacts = exact <$> ts
    let ys = evalSol de ts interps
    let max_err = maximum . fmap norm $ zipWith (^-^) ys exacts
    let last_err = norm $ last ys ^-^ last exacts
    return $
      if max_err / tol >= 1
        then testFailed ("unacceptable error " ++ show max_err ++ " | " ++ show last_err)
        else testPassed ""
  testOptions = pure []
