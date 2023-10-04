{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module Solver.ButcherTableau (BT (..), ERK (..), Tol (..), ERK_Params (..), dopri5, bosh3) where

import Control.Applicative (Const)
import Control.Monad.Reader (MonadReader (..))
import Control.Monad.State
import Data.Data (Proxy)
import Data.Maybe (fromJust)
import Data.Vector (Vector)
import qualified Data.Vector as V
import GHC.Exts (IsList)
import GHC.Generics (Generic)
import Interpolate
import Linear
import Linear.V
import Optics
import Solver.Class (ErrEst (..), Solver (..))
import Term

-- | Weighted sum of vector s
wsum :: (Num a, Additive v) => V.Vector a -> V.Vector (v a) -> v a
wsum cs vs = V.foldr (^+^) zero (V.zipWith (*^) cs vs)

data BT a = BT
  { coeffs :: ([(a, Vector a)]),
    f1 :: Vector a,
    f2 :: Vector a
  }
  deriving (Generic)

butcherTableau ::
  forall a ode.
  (Floating a, Term ode) =>
  [(a, Vector a)] ->
  ode a ->
  T ode a ->
  a ->
  a ->
  V.Vector (S ode a)
butcherTableau coeffs ode y0 t0 h = foldl go [] coeffs
  where
    go :: Vector (S ode a) -> (a, Vector a) -> Vector (S ode a)
    go ks (s, as) = snoc ks v
      where
        y = y0 ^+^ prod ode (wsum as ks) u
        t = t0 + s * h
        v = vf ode t y
        u = control ode (t0, t0 + h)

-- k = prod ode v u

-- | Explicit Runge Kutta
data ERK v a = ERK {bt :: BT a, interpCoeff :: Maybe (V.Vector a)}

data ERK_State v a = ERK_State {errEst :: a}
  deriving (Generic)

data ERK_Params v a = ERK_Params
  {tol :: Tol a}
  deriving (Generic)

data RKF45 v a = RKF45 (v a) (v a) -- Error estimate

data Tol a = Tol {aTol :: a, rTol :: a}

instance (Term ode) => Solver ERK ode where
  type SolState ERK ode = ERK_State (T ode)
  type SolParams ERK ode = ERK_Params (T ode)

  initSolver _ _ _ _ _ = ERK_State 0
  step (ERK (BT {..}) ics) ode y0 (t0, t1) = do
    Tol {..} <- view #tol <$> ask

    let ks = butcherTableau coeffs ode y0 t0 (t1 - t0)
    let u = control ode (t0, t1) -- is this the right thing??
    let y1 = y0 ^+^ prod ode (wsum f1 ks) u

    -- for error estimation
    let y1' = y0 ^+^ prod ode (wsum f2 ks) u
    let yy = liftI2 (\y y' -> max (abs y) (abs y')) y1 y1'
    let diff = y1 ^-^ y1'
    let tol = (rTol *^ yy) <&> (+ aTol)
    let err = norm $ liftI2 (/) diff tol
    assign #errEst err

    -- return $ H3 (t0, t1) (V.head ks) y0 (V.last ks) y1

    return $ case ics of
      Nothing -> H3 (t0, t1) (V.head ks) y0 (V.last ks) y1
      Just ics ->
        let ymid =
              y0
                ^+^ prod
                  ode
                  (wsum ics ks)
                  u -- (control ode (t0, t0 + 1 {- this or u? -}))
         in H4 (t0, t1) ymid y0 y1 (V.head ks) (V.last ks)

instance Term ode => ErrEst ERK ode where
  errorEstimate _ _ (ERK_State errEst) = errEst
  errorOrder de (ERK {bt}) = length $ f1 bt

--------

rkf45_bt :: (Floating a) => BT a
rkf45_bt =
  BT
    { f1 =
        -- fromJust . fromVector $
        [25 / 216, 0, 1408 / 2565, 2197 / 4104, -1 / 5, 0],
      f2 =
        -- fromJust . fromVector $
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

dopri_bt :: Floating a => BT a
dopri_bt =
  BT
    { f1 =
        -- fromJust . fromVector $
        [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0],
      f2 =
        -- fromJust . fromVector $
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

-- from page 6 of "Runge–Kutta pairs of orders 5(4) satisfying only the first column simplifying assumption"
-- http://users.uoa.gr/~tsitourasc/RK54_new_v2.pdf
tsit5_interp :: Floating a => Vector (Vector a)
tsit5_interp = undefined

dopri5_interp :: Floating a => Vector a
dopri5_interp =
  V.fromList
    [ 6025192743 / 30085553152 / 2,
      0,
      51252292925 / 65400821598 / 2,
      -2691868925 / 45128329728 / 2,
      187940372067 / 1594534317056 / 2,
      -1776094331 / 19743644256 / 2,
      11237099 / 235043384 / 2
    ]

{-# SPECIALIZE dopri5 :: ERK v Double #-}
dopri5 :: Floating a => ERK v a
dopri5 = ERK {bt = dopri_bt, interpCoeff = Just dopri5_interp}

bosh3_bt :: Floating a => BT a
bosh3_bt =
  BT
    { coeffs =
        [ (1 / 2, [1 / 2]),
          (3 / 4, [0.0, 3 / 4]),
          (1.0, [2 / 9, 1 / 3, 4 / 9])
        ],
      f1 = [2 / 9, 1 / 3, 4 / 9, 0.0],
      f2 = [7 / 24, 1 / 4, 1 / 3, -1 / 8]
    }

-- | Bogacki--Shampine's 3/2 method aka Ralston's third order
bosh3 :: Floating a => ERK v a
bosh3 = ERK bosh3_bt Nothing
