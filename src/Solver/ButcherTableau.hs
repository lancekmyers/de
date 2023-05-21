{-# LANGUAGE UndecidableInstances #-}

module Solver.ButcherTableau where

import Control.Monad.State
import Data.Data (Proxy)
import Linear
import Optics (snoc)
import Solver.Class
import Term

newtype BT a = BT ([a], [(a, [a])])

butcherTableau ::
  forall a ode.
  (Floating a, Term ode) =>
  [(a, [a])] ->
  ode a ->
  T ode a ->
  a ->
  a ->
  [T ode a]
butcherTableau coeffs ode y0 t0 h = foldl go [] coeffs
  where
    -- go :: (a, [a]) -> [T ode a] -> [T ode a]
    go ks (s, as) = snoc ks k
      where
        y = foldr (^+^) y0 (zipWith (*^) as ks)
        t = t0 + s * h
        v = vf ode t y
        u = control ode (t0, t0 + h)
        k = prod ode v u

rkf45 :: (Floating a) => ([a], [a], [(a, [a])])
rkf45 =
  ( [25 / 216, 0, 1408 / 2565, 2197 / 4104, -1 / 5, 0],
    [16 / 135, 0, 6656 / 12825, 28561 / 56430, -9 / 50, 2 / 55],
    [ (0, []),
      (0.25, [0.25]),
      (3 / 8, [3 / 32, 9 / 32]),
      (12 / 13, [1932 / 2197, -7200 / 2197, 7296 / 2197]),
      (1, [439 / 216, -8, 3680 / 513, -845 / 4104]),
      (0.5, [-8 / 27, 2, -3544 / 2565, 1859 / 4104, -11 / 40])
    ]
  )

data RKF45 v a = RKF45 (v a) (v a) -- Error estimate

instance
  (Term ode, Metric (T ode)) =>
  Solver RKF45 ode
  where
  initSolver _ _ _ = RKF45 zero zero
  step ode y0 (t0, t1) = do
    let (final1, final2, coeffs) = rkf45
    let ks = butcherTableau coeffs ode y0 t0 (t1 - t0)
    let y1 = foldr (^+^) y0 $ zipWith (*^) final1 ks
    let y1' = foldr (^+^) y0 $ zipWith (*^) final2 ks
    let yy = liftI2 (\y y' -> max (abs y) (abs y')) y1 y1'
    let err = y1 ^-^ y1'
    put $ RKF45 err yy
    return y1

instance
  (Term ode, Metric (T ode)) =>
  ErrEst RKF45 ode
  where
  errorOrder _ _ = 4
  errorEstimate atol rtol de (RKF45 err yy) =
    let tol = rtol *^ yy
     in norm $ liftI2 (\e t -> e / (atol + t)) err tol

stepButcher ::
  (Term de, Floating a, MonadState (sol (T de) a) m) =>
  BT a ->
  de a ->
  (T de) a ->
  (a, a) ->
  m (T de a)
stepButcher (BT (final, coeffs)) ode y0 (t0, t1) = do
  let ks = butcherTableau coeffs ode y0 t0 (t1 - t0)
  let y1 = foldr (^+^) y0 $ zipWith (*^) final ks
  return y1
