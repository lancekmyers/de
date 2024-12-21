module Interpolate
  ( Interp (..),
    interp,
    rightMost,
    evalSol,
    knots,
    mkH4,
    mkH3,
    mkLin,
  )
where

import Data.Vector (Vector)
import qualified Data.Vector as V
import Linear
import Optics (prism')
import Term
import Debug.Trace (traceShow)

data Interp v a
  = Poly (a, a) (Vector (v a))
  deriving Show

interp ::
  forall de a.
  (Term de, Floating a) =>
  de a ->
  a ->
  Interp (T de) a ->
  T de a
interp de c (Poly (t0, t1) coeffs) = polyEval coeffs c'
  where
    c' = (c - t0) / (t1 - t0)

mkLin ::
  (Additive v, Num a) =>
  de a ->
  (a, a) ->
  v a ->
  v a ->
  Interp v a
mkLin _de (t0, t1) y0 y1 = Poly (t0, t1) (V.fromList [y1 ^-^ y0, y0])

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
  Interp (T de) a
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

mkH3 ::
  (Term de, Floating a) =>
  de a ->
  (a, a) ->
  T de a ->
  S de a ->
  T de a ->
  S de a ->
  Interp (T de) a
mkH3 de (t0, t1) y0 f0 y1 f1 = Poly (t0, t1) coeffs
  where
    coeffs =
      V.fromList
        [ 2 *^ y0 ^+^ v0 ^-^ 2 *^ y1 ^+^ v1,
          (-3) *^ y0 ^-^ 2 *^ v0 ^+^ 3 *^ y1 ^-^ v1,
          v0,
          y0
        ]
    w = control de (t0, t1)
    v0 = prod de f0 w
    v1 = prod de f1 w

rightMost :: (Num a, Additive v) => Interp v a -> v a
rightMost (Poly _ coeffs) = V.foldl' (^+^) zero coeffs

leftMost :: (Num a, Additive v) => Interp v a -> v a
leftMost (Poly _ coeffs) = polyEval coeffs 0 -- V.last coeffs

timeInterval :: Interp v a -> (a, a)
timeInterval (Poly interval _) = interval

polyEval :: (Num a, Additive v) => V.Vector (v a) -> a -> v a
polyEval coeffs t = V.foldl' (\acc x -> x ^+^ (t *^ acc)) zero coeffs

contains :: (Ord a) => Interp v a -> a -> Bool
contains interp t =
  let (l, u) = timeInterval interp
   in t >= l && t <= u

type Solution de a = [Interp (T de) a]

-- | Evaluate solution at a list of times
-- requires times to be in ascending order
evalSol ::
  (Term de, Show a, Floating a, Ord a) =>
  de a ->
  [a] ->
  Solution de a ->
  [T de a]
evalSol de [] _ = []
evalSol de rest@(t:_) [] = error $ "not enough intervals " ++ show t ++ " " ++ show (length rest)
evalSol de (t : ts) (i : is)
  | contains i t = (interp de t i) : evalSol de ts (i : is)
  | otherwise = evalSol de (t : ts) is

knots :: (Num a, Additive v) => [Interp v a] -> [(a, v a)]
knots [] = []
knots (i : is) =
  (fst $ timeInterval i, leftMost i) :
    [(snd $ timeInterval i, rightMost i) | i <- is]
