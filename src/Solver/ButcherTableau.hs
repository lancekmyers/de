{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE UndecidableInstances #-}

module Solver.ButcherTableau where

import Control.Monad.State
import Data.Data (Proxy)
import Data.Maybe (fromJust)
import Data.Vector (Vector)
import qualified Data.Vector as V
import GHC.Exts (IsList)
import Linear
import Linear.V
import Optics (snoc)
import Solver.Class
import Term

newtype BT a = BT ([a], [(a, [a])])

data BT' n a = BT'
  { coeffs :: ([(a, Vector a)]),
    f1 :: V n a,
    f2 :: V n a
  }

butcherTableau ::
  forall a ode.
  (Floating a, Term ode) =>
  [(a, Vector a)] ->
  ode a ->
  T ode a ->
  a ->
  a ->
  V.Vector (T ode a)
butcherTableau coeffs ode y0 t0 h = foldl go [] coeffs
  where
    go :: Vector (T ode a) -> (a, Vector a) -> Vector (T ode a)
    go ks (s, as) = snoc ks k
      where
        y = foldr (^+^) y0 (V.zipWith (*^) as ks)
        t = t0 + s * h
        v = vf ode t y
        u = control ode (t0, t0 + h)
        k = prod ode v u

rkf45 :: (Floating a) => (BT' 6 a)
rkf45 =
  BT'
    { f1 =
        fromJust . fromVector $
          [25 / 216, 0, 1408 / 2565, 2197 / 4104, -1 / 5, 0],
      f2 =
        fromJust . fromVector $
          [16 / 135, 0, 6656 / 12825, 28561 / 56430, -9 / 50, 2 / 55],
      coeffs =
        [ (0, []),
          (0.25, [0.25]),
          (3 / 8, [3 / 32, 9 / 32]),
          (12 / 13, [1932 / 2197, -7200 / 2197, 7296 / 2197]),
          (1, [439 / 216, -8, 3680 / 513, -845 / 4104]),
          (0.5, [-8 / 27, 2, -3544 / 2565, 1859 / 4104, -11 / 40])
        ]
    }

data RKF45 v a = RKF45 (v a) (v a) -- Error estimate

instance
  (Term ode, Metric (T ode)) =>
  Solver RKF45 ode
  where
  initSolver _ _ _ = RKF45 zero zero
  step ode y0 (t0, t1) = do
    let BT' {..} = rkf45
    let ks = butcherTableau coeffs ode y0 t0 (t1 - t0)
    let y1 = foldr (^+^) y0 $ V.zipWith (*^) (toVector f1) ks
    let y1' = foldr (^+^) y0 $ V.zipWith (*^) (toVector f2) ks
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
  BT' n a ->
  de a ->
  (T de) a ->
  (a, a) ->
  m (T de a)
stepButcher BT' {..} ode y0 (t0, t1) = do
  let ks = butcherTableau coeffs ode y0 t0 (t1 - t0)
  let y1 = foldr (^+^) y0 $ V.zipWith (*^) (toVector f1) ks
  return y1

dopri :: BT' 7 Double
dopri =
  BT'
    { f1 =
        fromJust . fromVector $
          [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0],
      f2 =
        fromJust . fromVector $
          [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40],
      coeffs =
        [ (0, []),
          (1 / 5, [1 / 5]),
          (3 / 10, [3 / 40, 9 / 40]),
          (4 / 5, [44 / 45, -56 / 15, 32 / 9]),
          (8 / 9, [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729]),
          (1, [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656]),
          (1, [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84])
        ]
    }
