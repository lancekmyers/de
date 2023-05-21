module Main (main) where

import Linear
import Solver
import Term
import Test.Tasty (defaultMain, testGroup)
import Test.Tasty.Providers (singleTest)
import Tester

expIVP :: IVP (SimpleODE V1) Double
expIVP = IVP (SimpleODE (\t x -> x)) (\t -> V1 (exp t)) (0, 1) (V1 1)

testEuler ivp = DETest EulerSol (ConstantStep 5e-2) expIVP 1e-1

testHeun ivp = DETest HeunSol (ConstantStep 5e-2) expIVP 1e-2

testRKF ivp = DETest (RKF45 zero zero) (ConstantStep 1e-2) expIVP 1e-3

testRKFAdapt ivp = DETest (RKF45 zero zero) (SimpleAdaptStep 1e-4 1e-3) expIVP 1e-2

expTests =
  testGroup "exp" $
    [ singleTest "euler" $ testEuler expIVP,
      singleTest "heun" $ testHeun expIVP,
      singleTest "rkf45" $ testRKF expIVP,
      singleTest "rkf45 adapt" $ testRKFAdapt expIVP
    ]

main :: IO ()
main = defaultMain expTests
