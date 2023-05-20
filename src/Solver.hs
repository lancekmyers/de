{-# LANGUAGE DataKinds #-}
{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE UndecidableInstances #-}

module Solver (module BTab, module Euler, module C) where

import Solver.ButcherTableau as BTab
import Solver.Class as C
import Solver.Euler as Euler
