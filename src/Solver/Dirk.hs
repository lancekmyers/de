module Solver.Dirk where

import Control.Monad.Identity (Identity)
import Data.Foldable (Foldable (..))
import Data.Functor.Rep
import Data.Machine (auto, autoM)
import Data.Traversable (mapAccumL)
import Data.Vector qualified as V
import GHC.TypeNats (KnownNat)
import Interpolate
import Linear (Additive, zero, (*^), (^+^), (^-^))
import Linear.Affine (Vector)
import Linear.Matrix
import Linear.Metric (Metric)
import Linear.V
import Numeric.AD.Rank1.Forward (Forward, jacobian)
import Solver.Class (ErrorEstimate, Solver, TimeStep (..), Tol (..))
import Term

wsum :: (Num a, Additive v) => V.Vector a -> V.Vector (v a) -> v a
wsum cs vs = V.foldr (^+^) zero (V.zipWith (*^) cs vs)

data BT n a = BT
  { c :: V n a,
    a :: V n (V.Vector a),
    b :: V n a,
    b_hat :: V n a,
    interp :: IC a
  }

-- erk :: forall ode a. (Floating a, Ord a, Term ode) => BT a -> IC a -> Tol a -> ode a -> Solver (ErrorEstimate a) Identity (T ode) a
dirk ::
  forall ode v a n.
  ( Floating a,
    Ord a,
    Metric v,
    Finite v,
    KnownNat (Size v),
    Num (v a),
    Traversable v,
    Applicative v,
    Foldable v
  ) =>
  BT n a ->
  JacODE v a ->
  Solver () Identity (T (JacODE v)) a
dirk bt ode = autoM $ pure . go
  where
    go (y0, tstp) =
      let y1 = step bt ode y0 tstp
       in let TimeStep {t, delta} = tstp
           in ((), mkLin ode (t, t + delta) y0 y1)

step ::
  forall v a n.
  ( Floating a,
    Finite v,
    Metric v,
    Foldable v,
    Traversable v,
    Applicative v,
    KnownNat (Size v),
    Num (v a)
  ) =>
  BT n a ->
  JacODE v a ->
  v a ->
  TimeStep a ->
  v a
step bt ode y0 (TimeStep t0 h) = wsum (toVector b) ks
  where
    j = jac ode t0 y0
    BT {c, a, b, b_hat} = bt
    ts = V.map (\c -> t0 + h * c) (toVector c)
    go :: V.Vector (v a) -> (a, V.Vector a) -> V.Vector (v a)
    go acc (ti, ai) = V.snoc acc $ (stage (vf ode ti) j xi (h * a_diag) y0) !! 2
      where
        Just (a_lower, a_diag) = V.unsnoc ai
        xi :: v a
        xi = wsum a_lower acc
    ks = foldl' go V.empty (V.zip ts (toVector a))

stage ::
  forall v a.
  ( Fractional a,
    Finite v,
    Additive v,
    Foldable v,
    Traversable v,
    Applicative v,
    KnownNat (Size v),
    Num (v a)
  ) =>
  (v a -> v a) ->
  v (v a) ->
  v a ->
  a ->
  v a ->
  [v a]
stage f jac x a_diag y = iterate go y
  where
    r y' = x ^+^ a_diag *^ (f y') ^-^ (y' ^-^ y)
    newt_iter :: v (v a)
    newt_iter = identity !-! a_diag *!! jac
    newt_iter_inv = luInvFinite newt_iter
    go :: v a -> v a
    go u = u ^+^ newt_iter_inv !* (r u)
