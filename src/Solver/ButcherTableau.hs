module Solver.ButcherTableau where

import Control.Monad.State
import Linear
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
butcherTableau coeffs ode y0 t0 h = foldr go [] coeffs
  where
    go :: (a, [a]) -> [T ode a] -> [T ode a]
    go (s, as) ks = k : ks
      where
        y = foldr (^+^) zero (zipWith (*^) as ks)
        t = t0 + s * h
        v = vf ode t y
        u = control ode (t0, t0 + h)
        k = prod ode v u

rk45 :: Floating a => BT a
rk45 =
  BT
    ( [1 / 6, 1 / 3, 1 / 3, 1 / 6],
      [ (0, []),
        (0.5, [0.5]),
        (0.5, [0, 0.5]),
        (1.0, [0, 0, 1])
      ]
    )

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
