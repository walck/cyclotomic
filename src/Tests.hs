{-# OPTIONS_GHC -Wall #-}

module Main where

import Data.Complex.Cyclotomic
import Test.Framework
    ( defaultMain
    , testGroup
    )
import Test.Framework.Providers.HUnit
    ( testCase
    )
import Test.Framework.Providers.QuickCheck2
    ( testProperty
    )
import qualified Test.Framework.Providers.SmallCheck as S
import qualified Test.Framework.Providers.API as T
import Test.QuickCheck
    ( Gen
    , elements
    , Arbitrary(..)
    , shrinkRealFrac
    )
import Test.HUnit
    ( (@?=)
    , Assertion
    )
import Data.List
    ( nub
    )
import Data.Ratio
    ( (%)
    )

main :: IO ()
main = defaultMain tests

tests :: [T.Test]
tests = [test1a
        ,test2b
        ,test3b
        ,test4b
        ,S.withDepth 10 (S.testProperty "SmallCheck prop_square_sqrtRat" prop_square_sqrtRat)
        ,qc_square_sqrtRat
        ,qc_Gauss
        ,qc_dftInv_dft
        ,qc_dft_dftInv
        ,qc_sum_quadratic_roots
        ]

rationals :: [Rational]
rationals = 0 % 1 : [sign * k % j | n <- [0..], m <- [0..n-1], sign <- [1,-1]
                    , let k = m + 1, let j = n - m, gcd k j == 1]

rationalList :: Integer -> [Rational]
rationalList m = nub [n % d | n <- [-m..m], d <- [1..m]]

test1a :: T.Test
test1a = testGroup "polarRat" [polarRatTest p q | p <- [0..10], q <- [1..10]]

polarRatAssertion :: Integer -> Integer -> Assertion
polarRatAssertion p q = polarRat 1 (p % q) @?= e q^p

polarRatTest :: Integer -> Integer -> T.Test
polarRatTest p q = testCase ("polarRat 1 (" ++ show p ++ " % " ++ show q ++ ")") (polarRatAssertion p q)

test2b :: T.Test
test2b = testGroup "sqrtRat r ^ 2 == r for the following values of r"
         [testCase (show r) (sqrtRat r ^ (2::Int) @?= fromRational r)
              | r <- take 100 rationals]

test3b :: T.Test
test3b = testGroup "sqrtRat (r*r) == abs r for the following values of r"
         [testCase (show r) (sqrtRat (r*r) @?= fromRational (abs r))
              | r <- take 100 rationals]

test4b :: T.Test
test4b = testGroup "z * (1 / z) == 1 for the following values of z"
         [testCase (show z) (z * (1 / z) @?= 1)
              | n <- [1..10], m <- [1..10], let z = e n + e m, z /= 0]

----------------
-- Properties --
----------------

prop_square_sqrtRat :: Int -> Bool
prop_square_sqrtRat n = sqrtRat (fromIntegral n) ^ (2::Int) == fromIntegral n

prop_Gauss :: Integer -> Bool
prop_Gauss n = let nn = 2 * abs n + 1
               in sum [e nn^(j*j `mod` nn) | j <- [1..(nn - 1) `div` 2]]
                      == if nn `mod` 4 == 1
                         then (-1 + sqrtInteger nn) / 2
                         else (-1 + i*sqrtInteger nn) / 2

prop_dftInv_dft :: [Rational] -> Bool
prop_dftInv_dft rs = dftInv (dft cs) == cs
    where cs = map fromRational rs

prop_dft_dftInv :: [Rational] -> Bool
prop_dft_dftInv rs = dft (dftInv cs) == cs
    where cs = map fromRational rs

prop_sum_quadratic_roots :: (Rational, Rational, Rational) -> Bool
prop_sum_quadratic_roots (a, b, c)
    = case rootsQuadEq a b c of
        Nothing      -> a == 0
        Just (r1,r2) -> r1 + r2 == fromRational (-b / a)

prop_sum_quadratic_roots_small :: (SmallRational, SmallRational, SmallRational) -> Bool
prop_sum_quadratic_roots_small (SmallRational a, SmallRational b, SmallRational c)
    = case rootsQuadEq a b c of
        Nothing      -> a == 0
        Just (r1,r2) -> r1 + r2 == fromRational (-b / a)

----------------------
-- QuickCheck Tests --
----------------------

qc_square_sqrtRat :: T.Test
qc_square_sqrtRat
    = T.plusTestOptions (T.TestOptions
                         {T.topt_seed                               = Nothing
                         ,T.topt_maximum_generated_tests            = Just 15
                         ,T.topt_maximum_unsuitable_generated_tests = Nothing
                         ,T.topt_maximum_test_size                  = Just 15
                         ,T.topt_maximum_test_depth                 = Nothing
                         ,T.topt_timeout                            = Nothing
                         })
      $ testProperty "QuickCheck prop_square_sqrtRat" prop_square_sqrtRat

qc_Gauss :: T.Test
qc_Gauss
    = T.plusTestOptions (T.TestOptions
                         {T.topt_seed                               = Nothing
                         ,T.topt_maximum_generated_tests            = Nothing
                         ,T.topt_maximum_unsuitable_generated_tests = Nothing
                         ,T.topt_maximum_test_size                  = Nothing
                         ,T.topt_maximum_test_depth                 = Nothing
                         ,T.topt_timeout                            = Nothing
                         })
      $ testProperty "QuickCheck prop_Gauss" prop_Gauss

qc_dftInv_dft :: T.Test
qc_dftInv_dft
    = T.plusTestOptions (T.TestOptions
                         {T.topt_seed                               = Nothing
                         ,T.topt_maximum_generated_tests            = Just 15
                         ,T.topt_maximum_unsuitable_generated_tests = Nothing
                         ,T.topt_maximum_test_size                  = Just 30
                         ,T.topt_maximum_test_depth                 = Nothing
                         ,T.topt_timeout                            = Nothing
                         })
      $ testProperty "QuickCheck prop_dftInv_dft" prop_dftInv_dft

qc_dft_dftInv :: T.Test
qc_dft_dftInv
    = T.plusTestOptions (T.TestOptions
                         {T.topt_seed                               = Nothing
                         ,T.topt_maximum_generated_tests            = Just 15
                         ,T.topt_maximum_unsuitable_generated_tests = Nothing
                         ,T.topt_maximum_test_size                  = Just 30
                         ,T.topt_maximum_test_depth                 = Nothing
                         ,T.topt_timeout                            = Nothing
                         })
      $ testProperty "QuickCheck prop_dft_dftInv" prop_dft_dftInv

qc_sum_quadratic_roots :: T.Test
qc_sum_quadratic_roots
    = testProperty "QuickCheck prop_sum_quadratic_roots" prop_sum_quadratic_roots_small

----------------------
-- QuickCheck Stuff --
----------------------

data SmallRational = SmallRational Rational
                     deriving (Show,Ord,Eq)

smallRationalList :: [SmallRational]
smallRationalList = map SmallRational (rationalList 3)

smallRationalGen :: Gen SmallRational
smallRationalGen = elements smallRationalList

instance Arbitrary SmallRational where
    arbitrary = smallRationalGen
    shrink (SmallRational r) = map SmallRational (shrinkRealFrac r)

