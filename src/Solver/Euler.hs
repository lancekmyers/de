{-# LANGUAGE TypeFamilies #-}

module Solver.Euler (euler) where

import Control.Applicative (Const (..))
import Control.Monad.Reader
import Control.Monad.State
import Data.Machine.Mealy
import Interpolate
import Linear
import Solver.Class
import Term

euler ::
  (Ord a, Floating a, Term ode) =>
  ode a ->
  StepIntegrator (T ode) a
euler ode = unfoldMealy go ()
  where
    go _ (y, TimeStep {t = t0, delta}) =
      let t1 = t0 + delta
          y' = vf ode t0 y
          dt = control ode (t0, t1)
          y1 = y ^+^ prod ode y' dt
          interp = mkLin ode (t0, t1) y y1
       in ((interp, undefined), ())

-- data HeunSol v a = HeunSol

-- instance (Term ode) => Solver HeunSol ode where
--   type SolState HeunSol ode = Const ()
--   type SolParams HeunSol ode = Const ()

--   initSolver _ _ _ _ _ = Const ()
--   step _sol ode y (t0, t1) = do
--     -- EulerStep
--     let y1 =
--           rightMost $
--             step EulerSol ode y (t0, t1)
--               `evalStateT` Const ()
--               `runReader` Const ()
--     let v1 = vf ode t0 y
--     let v2 = vf ode t1 y1
--     -- average the velocity
--     let v = (v1 ^+^ v2) ^/ 2
--     let dt = control ode (t0, t1)
--     let y2 = y ^+^ prod ode v dt
--     return $ mkLin ode (t0, t1) y y2
