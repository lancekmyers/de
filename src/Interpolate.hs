module Interpolate (Interp (..), interp, rightMost, evalSol, knots, mkH4) where

import Data.Vector (Vector)
import qualified Data.Vector as V
import Linear
import Optics (prism')
import Term

data Interp de a
  = H3
      { interval :: (a, a),
        m :: S de a,
        p :: T de a,
        m' :: S de a,
        p' :: T de a
      }
  | Lin (a, a) (T de a) (T de a)
  | Poly (a, a) (Vector (T de a))

interp ::
  forall de a.
  (Term de, Floating a) =>
  de a ->
  a ->
  Interp de a ->
  T de a
interp de c (Poly (t0, t1) coeffs) = polyEval coeffs c'
  where
    c' = (c - t0) / (t1 - t0)
interp de x (Lin (x0, x1) p p') =
  t *^ p ^+^ (1 - t) *^ p'
  where
    t = (x - x0) / x1 - x0
interp de x (H3 (x0, x1) m p m' p') =
  h00 *^ p
    ^+^ h10 *^ (prod de m w)
    ^+^ h01 *^ p'
    ^+^ h11 *^ (prod de m' w)
  where
    w :: U de a
    w = control de (x0, x1)
    t = (x - x0) / x1 - x0
    h00 = (1 + 2 * t) * (1 - t) ^^ 2
    h10 = t * (1 - t) ^^ 2
    h01 = t ^^ 2 * (3 - 2 * t)
    h11 = t ^^ 2 * (t - 1)

-- | Fourth Order Hermite polynomial interpolation
mkH4 ::
  forall de a.
  (Term de, Floating a) =>
  de a ->
  (a, a) ->
  T de a ->
  T de a ->
  T de a ->
  S de a ->
  S de a ->
  Interp de a
mkH4 de (t0, t1) ymid y0 y1 f0 f1 = Poly (t0, t1) coeffs
  where
    w = control de (t0, t1) -- control de (0, 1)
    f1' = prod de f1 w
    f0' = prod de f0 w
    _a =
      2 *^ (f1' ^-^ f0')
        ^-^ 8 *^ (y1 ^+^ y0) ^+^ 16 *^ ymid
    _b =
      5 *^ f0' ^-^ 3 *^ f1'
        ^+^ 18 *^ y0
        ^+^ 14 *^ y1 ^-^ 32 *^ ymid
    _c =
      f1' ^-^ 4 *^ f0'
        ^-^ 11 *^ y0
        ^-^ 5 *^ y1 ^+^ 16 *^ ymid
    coeffs = V.fromList [_a, _b, _c, f0', y0]

rightMost :: (Num a, Additive (T de)) => Interp de a -> T de a
rightMost (Poly _ coeffs) = V.foldl' (^+^) zero coeffs
rightMost (Lin _ _ y1) = y1
rightMost (H3 _ _ _ _ y1) = y1

leftMost :: (Num a, Additive (T de)) => Interp de a -> T de a
leftMost (Poly _ coeffs) = polyEval coeffs 0 -- V.last coeffs
leftMost (Lin _ y0 _) = y0
leftMost (H3 _ _ y0 _ _) = y0

timeInterval :: Interp de a -> (a, a)
timeInterval (Poly interval _) = interval
timeInterval (H3 interval _ _ _ _) = interval
timeInterval (Lin interval _ _) = interval

polyEval :: (Num a, Additive v) => V.Vector (v a) -> a -> v a
polyEval coeffs t = V.foldl' (\acc x -> x ^+^ (t *^ acc)) zero coeffs

contains :: (Ord a) => Interp de a -> a -> Bool
contains interp t =
  let (l, u) = timeInterval interp
   in t >= l && t <= u

type Solution de a = [Interp de a]

-- | Evaluate solution at a list of times
-- requires times to be in ascending order
evalSol ::
  (Term de, Floating a, Ord a) =>
  de a ->
  [a] ->
  Solution de a ->
  [T de a]
evalSol de [] _ = []
evalSol de _ [] = []
evalSol de (t : ts) (i : is)
  | contains i t = interp de t i : evalSol de ts is
  | otherwise = evalSol de (t : ts) is

knots :: (Num a, Term de) => [Interp de a] -> [(a, T de a)]
knots [] = []
knots (i : is) =
  (fst $ timeInterval i, leftMost i) :
    [(snd $ timeInterval i, rightMost i) | i <- is]
