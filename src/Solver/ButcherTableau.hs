{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedLabels #-}
{-# LANGUAGE OverloadedLists #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module Solver.ButcherTableau (BT (..), ERK (..), Tol (..), ERK_Params (..), dopri5, bosh3, tsit5) where

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
data ERK v a = ERK {bt :: BT a, interpCoeff :: IC a}

-- | Iterpolation Coefficients
data IC a
  = InterpLinear
  | -- use the normal Hermite 3rd order polynomial
    PlainH3
  | -- coeffs for computing midpoint to construct Hermite 4th order polynomial
    MidPointH4 (V.Vector a)
  | InterpMatrix (V.Vector (V.Vector a))

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

    return $ case ics of
      InterpLinear -> mkLin ode (t0, t1) y0 y1
      PlainH3 -> mkH3 ode (t0, t1) y0 (V.head ks) y1 (V.last ks)
      MidPointH4 ics ->
        let ymid =
              y0
                ^+^ prod
                  ode
                  (wsum ics ks)
                  u
         in mkH4 ode (t0, t1) ymid y0 y1 (V.head ks) (V.last ks)
      InterpMatrix mat ->
        let coeffs = ((flip (prod ode) u) . (flip wsum ks) <$> mat) `V.snoc` y0
         in Poly (t0, t1) coeffs

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
dopri5 = ERK {bt = dopri_bt, interpCoeff = MidPointH4 dopri5_interp}

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
bosh3 = ERK bosh3_bt PlainH3

-- from page 6 of "Runge–Kutta pairs of orders 5(4) satisfying only the first column simplifying assumption"
-- http://users.uoa.gr/~tsitourasc/RK54_new_v2.pdf
tsit5_interp :: Floating a => Vector (Vector a)
tsit5_interp =
  [ [ -1.053088497729022,
      0.101700000000000,
      2.49062728565125279,
      -16.5481028892449027,
      47.3795219628192812,
      -34.8706578614966097,
      2.50000000000000
    ], -- t^4
    [ 2.913255461821913,
      -0.223400000000000,
      -5.9410338721315047,
      30.33818863028232,
      -88.178904894766401,
      65.0918946747937,
      -4.00000000000000
    ], -- t^3
    [ -2.76370619727483,
      0.131700000000000,
      3.9302962368947515,
      -12.41107716693368,
      37.509313416511039,
      -27.8965262891973,
      1.50000000000000
    ], -- t^2
    [1, 0, 0, 0, 0, 0, 0], -- t^1
    zero -- t^0
  ]

-- https://github.com/patrick-kidger/diffrax/blob/9126ce0c7656951945a8b779722e19b7ebddca33/diffrax/solver/tsit5.py#L38C14-L38C14
tsit5_bt :: Floating a => BT a
tsit5_bt =
  BT
    { coeffs =
        [ (c0,) [161 / 1000],
          (c1,)
            [ -0.8480655492356988544426874250230774675121177393430391537369234245294192976164141156943e-2,
              0.3354806554923569885444268742502307746751211773934303915373692342452941929761641411569
            ],
          (c2,)
            [ 2.897153057105493432130432594192938764924887287701866490314866693455023795137503079289,
              -6.359448489975074843148159912383825625952700647415626703305928850207288721235210244366,
              4.362295432869581411017727318190886861027813359713760212991062156752264926097707165077
            ],
          (c3,)
            [ 5.325864828439256604428877920840511317836476253097040101202360397727981648835607691791,
              -11.74888356406282787774717033978577296188744178259862899288666928009020615663593781589,
              7.495539342889836208304604784564358155658679161518186721010132816213648793440552049753,
              -0.9249506636175524925650207933207191611349983406029535244034750452930469056411389539635e-1
            ],
          (c4,)
            [ 5.861455442946420028659251486982647890394337666164814434818157239052507339770711679748,
              -12.92096931784710929170611868178335939541780751955743459166312250439928519268343184452,
              8.159367898576158643180400794539253485181918321135053305748355423955009222648673734986,
              -0.7158497328140099722453054252582973869127213147363544882721139659546372402303777878835e-1,
              -0.2826905039406838290900305721271224146717633626879770007617876201276764571291579142206e-1
            ],
          (c5,)
            [ 0.9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1,
              1 / 100,
              0.4798896504144995747752495322905965199130404621990332488332634944254542060153074523509,
              1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331,
              -3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677,
              2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841
            ]
        ],
      f1 =
        [ 0.9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1,
          1 / 100,
          0.4798896504144995747752495322905965199130404621990332488332634944254542060153074523509,
          1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331,
          -3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677,
          2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841,
          0.0
        ],
      f2 =
        [ 0.9646076681806522951816731316512876333711995238157997181903319145764851595234062815396e-1
            - 0.9468075576583945807478876255758922856117527357724631226139574065785592789071067303271e-1,
          1 / 100
            - 0.9183565540343253096776363936645313759813746240984095238905939532922955247253608687270e-2,
          0.4798896504144995747752495322905965199130404621990332488332634944254542060153074523509
            - 0.4877705284247615707855642599631228241516691959761363774365216240304071651579571959813,
          1.379008574103741893192274821856872770756462643091360525934940067397245698027561293331
            - 1.234297566930478985655109673884237654035539930748192848315425833500484878378061439761,
          -3.290069515436080679901047585711363850115683290894936158531296799594813811049925401677
            + 2.707712349983525454881109975059321670689605166938197378763992255714444407154902012702,
          2.324710524099773982415355918398765796109060233222962411944060046314465391054716027841
            - 1.866628418170587035753719399566211498666255505244122593996591602841258328965767580089,
          -1 / 66
        ]
    }
  where
    c0 = 161 / 1000
    c1 = 327 / 1000
    c2 = 9 / 10
    c3 = 0.9800255409045096857298102862870245954942137979563024768854764293221195950761080302604
    c4 = 1.0
    c5 = 1.0

tsit5 :: ERK v Double
tsit5 = ERK tsit5_bt (InterpMatrix tsit5_interp)
