module Main (main) where

import Control.Applicative (Const (..))
import Linear
import Solver
import Term
import Test.Tasty (defaultMain, testGroup)
import Test.Tasty.Providers (singleTest)
import Tester

otherIVP =
  IVP
    (SimpleODE $ \t y -> sin (t) *^ y)
    (\t -> V1 $ exp (1 - cos t))
    (0, 5)
    [x / 10 | x <- [0 .. 50]]
    (V1 1)

expIVP :: IVP (SimpleODE V1) Double
expIVP =
  IVP
    (SimpleODE (\t x -> x))
    (\t -> V1 (exp t))
    (0, 1)
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    (V1 1)

testEuler ivp = DETest EulerSol (Const ()) ConstantStepper ivp 1e-1

testHeun ivp = DETest HeunSol (Const ()) ConstantStepper ivp 1e-2

testDopri :: (Floating a, Real a, Show a, Ord a, Term ode) => IVP ode a -> DETest
testDopri ivp = DETest dopri5 (ERK_Params (Tol 1e-6 1e-6)) pi42 ivp 1e-3

testDopriAdapt ivp =
  DETest
    dopri5
    (ERK_Params (Tol 1e-3 1e-6))
    pi34
    ivp
    1e-4

testTsitAdapt :: Term ode => IVP ode Double -> DETest
testTsitAdapt ivp =
  DETest
    tsit5
    (ERK_Params (Tol {rTol = 1e-3, aTol = 1e-6}))
    h312PID
    ivp
    1e-4

testBoshAdapt :: Term ode => IVP ode Double -> DETest
testBoshAdapt ivp =
  DETest
    bosh3
    (ERK_Params (Tol 1e-4 1e-4))
    basicI
    ivp
    1e-2

expTests =
  testGroup "exp" $
    [ singleTest "euler" $ testEuler otherIVP,
      singleTest "heun" $ testHeun otherIVP,
      singleTest "dopri5-fixed" $ testDopri otherIVP,
      singleTest "dopri5-adapt" $ testDopriAdapt otherIVP,
      singleTest "bosh3-adapt" $ testBoshAdapt otherIVP,
      singleTest "tsit5-adapt" $ testTsitAdapt otherIVP
    ]

main :: IO ()
main = defaultMain expTests
