module Solver.Euler (EulerSol, HeunSol) where

import Control.Monad.State
import Linear
import Solver.Class
import Term

data EulerSol v a = EulerSol

instance (Term ode) => Solver EulerSol ode where
  initSolver _ _ _ = EulerSol
  step ode y (t0, t1) = do
    let y' = vf ode t0 y
    let dt = control ode (t0, t1)
    return $ y ^+^ prod ode y' dt

data HeunSol v a = HeunSol

instance (Term ode) => Solver HeunSol ode where
  initSolver _ _ _ = HeunSol
  step ode y (t0, t1) = do
    -- EulerStep
    let y1 = evalState (step ode y (t0, t1)) EulerSol
    let v1 = vf ode t0 y
    let v2 = vf ode t1 y1
    -- average the velocity
    let v = (v1 ^+^ v2) ^/ 2
    let dt = control ode (t0, t1)
    return $ y ^+^ prod ode v dt
