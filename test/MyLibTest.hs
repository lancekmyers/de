module Main (main) where

import Control.Applicative (Const (..))
import Linear
import Solver
import Term
import Test.Tasty (defaultMain, testGroup)
import Test.Tasty.Providers (singleTest)
import Tester

expIVP :: IVP (SimpleODE V1) Double
expIVP =
  IVP
    (SimpleODE (\t x -> x))
    (\t -> V1 (exp t))
    (0, 1)
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    (V1 1)

testEuler ivp = DETest EulerSol (Const ()) (ConstantStep 1e-2) expIVP 1e-1

testHeun ivp = DETest HeunSol (Const ()) (ConstantStep 1e-3) expIVP 1e-2

testDopri ivp = DETest dopri5 (ERK_Params (Tol 1e-6 1e-6)) (ConstantStep 1e-2) expIVP 1e-3

testDopriAdapt ivp =
  DETest
    dopri5
    (ERK_Params (Tol 1e-6 1e-6))
    SimpleAdaptStep
    expIVP
    1e-2

testBoshAdapt ivp =
  DETest
    bosh3
    (ERK_Params (Tol 1e-6 1e-6))
    SimpleAdaptStep
    expIVP
    1e-2

expTests =
  testGroup "exp" $
    [ singleTest "euler" $ testEuler expIVP,
      singleTest "heun" $ testHeun expIVP,
      singleTest "dopri5-fixed" $ testDopri expIVP,
      singleTest "dopri5-adapt" $ testDopriAdapt expIVP,
      singleTest "bosh3-adapt" $ testBoshAdapt expIVP
    ]

main :: IO ()
main = defaultMain expTests
