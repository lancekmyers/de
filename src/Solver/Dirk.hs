module Solver.Dirk where

import Data.Functor.Rep
import Data.Vector qualified as V
import GHC.TypeNats (KnownNat)
import Linear (Additive, zero, (*^), (^+^), (^-^))
import Linear.Matrix
import Linear.V
import Numeric.AD.Rank1.Forward (Forward, jacobian)

wsum :: (Num a, Additive v) => V.Vector a -> V.Vector (v a) -> v a
wsum cs vs = V.foldr (^+^) zero (V.zipWith (*^) cs vs)

data BT n a = BT
  { c :: V n a,
    a :: V n (V.Vector a),
    b :: V n a,
    b_hat :: V n a
  }


step ::
  forall v a n.
  ( Fractional a,
    Finite v,
    Additive v,
    Foldable v,
    Traversable v,
    Applicative v,
    KnownNat (Size v),
    Num (v a)
  ) =>
  BT n a ->
  JacODE v a ->
  v a ->
  (a, a) ->
  v a
step bt ode y0 (t0, h) = undefined
  where
    j = jac ode t0 y0
    BT {c, a, b} = bt
    ts = V.map (\c -> t0 + h * c) (toVector c)
    go acc (ti, ai) = V.snoc acc $ (stage (vf ode ti) j xi (h * a_diag) y0) !! 2
      where
        Just (a_lower, a_diag) = V.unsnoc ai
        xi :: v a
        xi = wsum a_lower acc

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
