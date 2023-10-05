{-# LANGUAGE DataKinds #-}
{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module Solver (module BTab, module Euler, module C) where

import Solver.ButcherTableau as BTab
  ( BT (..),
    ERK (..),
    ERK_Params (..),
    Tol (..),
    bosh3,
    dopri5,
    tsit5,
  )
import Solver.Class as C
  ( ConstantStep (..),
    ErrEst (..),
    SimpleAdaptStep (..),
    Solver (..),
    Stepper (..),
    solve,
  )
import Solver.Euler as Euler
