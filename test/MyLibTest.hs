module Main (main) where

import Control.Applicative (Const (..))
import Data.Typeable (Typeable)
import Linear
import Solver
import Term
import Test.Tasty (TestTree, defaultMain, testGroup)
import Test.Tasty.Providers (singleTest)
import Tester

otherIVP =
  IVP
    (SimpleODE $ \t y -> sin t *^ y)
    (\t -> V1 $ exp (1 - cos t))
    (0, 5)
    [x / 10 | x <- [0 .. 49]]
    (V1 1)

expIVP :: IVP (SimpleODE V1) Double
expIVP =
  IVP
    (SimpleODE (\t x -> x))
    (V1 . exp)
    (0, 1)
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    (V1 1)

testEuler ivp = DETest euler ivp 5e-2 1e-1

-- testHeun ivp = DETest heun (Const ()) ConstantStepper ivp 1e-2

testDopri :: (Typeable ode, Term ode) => IVP ode Double -> DETest ode Double
testDopri ivp = DETest (controlledStep pi33 . dopri5 (Tol 1e-4 1e-6)) ivp 1e-1 1e-3

testDopriAdapt ivp =
  DETest
    (controlledStep basicI . dopri5 (Tol 1e-5 1e-7))
    ivp
    1e-1
    1e-5

testTsitAdapt :: (Typeable ode, Term ode) => IVP ode Double -> DETest ode Double
testTsitAdapt ivp =
  DETest
    (controlledStep pi42 . tsit5 (Tol {rTol = 1e-4, aTol = 1e-6}))
    ivp
    1e-1
    1e-3

-- testBoshAdapt :: (Term ode) => IVP ode Double -> DETest
testBoshAdapt ivp =
  DETest
    (controlledStep basicI . bosh3 (Tol 1e-4 1e-6))
    ivp
    1e-1
    1e-3

testRKFAdapt ivp =
  DETest
    (controlledStep basicI . rkf45 (Tol 1e-4 1e-6))
    ivp
    1e-1
    5e-2

testIVP :: (Typeable ode, Term ode) => String -> IVP ode Double -> TestTree
testIVP name ivp =
  testGroup
    name
    [ singleTest "dopri5-fixed" $ testDopri ivp,
      singleTest "dopri5-adapt" $ testDopriAdapt ivp,
      -- bosh3 has borked error estimate
      singleTest "bosh3-adapt" $ testBoshAdapt ivp,
      singleTest "rkf45-adapt" $ testRKFAdapt ivp,
      singleTest "tsit5-adapt" $ testTsitAdapt ivp
    ]

main :: IO ()
main = defaultMain (testIVP "exp" expIVP)
