{-# LANGUAGE DataKinds #-}
{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module Solver (module BTab, module Euler, module C) where

import Solver.ButcherTableau as BTab
  ( BT (..),
    Tol (..),
    bosh3,
    dopri5,
    tsit5,
    rkf45
  )
import Solver.Class as C
  ( ErrorEstimate (..),
    TimeStep (..),
    basicI,
    constantStepper,
    h211PI,
    h312PID,
    pi33,
    pi34,
    pi42,
    runIntegration,
    solvingMachine,
    type StepController,
    type StepIntegrator,
  )
import Solver.Euler as Euler
