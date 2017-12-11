{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE Safe #-}

{- |
Module      :  Data.Number.RealCyclotomic
Copyright   :  (c) Scott N. Walck 2012-2017
License     :  GPL-3 (see LICENSE)
Maintainer  :  Scott N. Walck <walck@lvc.edu>
Stability   :  experimental

The real cyclotomic numbers are a subset of the real numbers with
the following properties:

     1.  The real cyclotomic numbers are represented exactly, enabling exact
     computations and equality comparisons.

     2.  The real cyclotomic numbers contain the rationals.
     As a consequence, the real cyclotomic numbers are a dense subset of the
     real numbers.

     3.  The real cyclotomic numbers contain the square roots of all nonnegative rational numbers.

     4.  The real cyclotomic numbers form a field:  they are closed under addition, subtraction,
     multiplication, and division.

     5.  The real cyclotomic numbers contain the sine and cosine of all rational
     multiples of pi (equivalently, the sine and cosine of any rational number
     of degrees or any rational number of revolutions).

     Floating point numbers do not do well with equality comparison:

>(sqrt 2 + sqrt 3)^2 == 5 + 2 * sqrt 6
> -> False

     "Data.Number.RealCyclotomic" represents these numbers exactly, allowing equality comparison:

>(sqrtRat 2 + sqrtRat 3)^2 == 5 + 2 * sqrtRat 6
> -> True

     'RealCyclotomic's can be exported as inexact real numbers using the 'toReal' function:

>sqrtRat 2
> -> e(8) - e(8)^3
>toReal $ sqrtRat 2
> -> 1.414213562373095

This module is based on the module 'Data.Complex.Cyclotomic'.
Usually you would only import one of the modules 'Data.Number.RealCyclotomic'
or 'Data.Complex.Cyclotomic', depending on whether you wanted only
real numbers (this module) or complex numbers (the other).
Functions such as @sqrtRat@, @sinDeg@, @cosDeg@ are defined
in both modules, with different type signatures, so their
names will conflict if both modules are imported.
-}

module Data.Number.RealCyclotomic
    ( RealCyclotomic
    , sqrtRat
    , sinDeg
    , cosDeg
    , sinRev
    , cosRev
    , isRat
    , toRat
    , toReal
    , goldenRatio
    , heron
    )
    where

import qualified Data.Complex.Cyclotomic as Cyc
import Data.Complex
    ( realPart
    )

-- | A real cyclotomic number.
newtype RealCyclotomic = RealCyclotomic Cyc.Cyclotomic
    deriving (Eq)

-- | @abs@ and @signum@ are undefined.
instance Num RealCyclotomic where
    RealCyclotomic x + RealCyclotomic y = RealCyclotomic (x + y)
    RealCyclotomic x - RealCyclotomic y = RealCyclotomic (x - y)
    RealCyclotomic x * RealCyclotomic y = RealCyclotomic (x * y)
    negate (RealCyclotomic x) = RealCyclotomic (negate x)
    fromInteger n = RealCyclotomic (fromInteger n)
    abs    = undefined
    signum = undefined

instance Fractional RealCyclotomic where
    recip (RealCyclotomic x) = RealCyclotomic (recip x)
    fromRational r = RealCyclotomic (fromRational r)

instance Show RealCyclotomic where
    show (RealCyclotomic x) = show x

-- I need to do Ord first.
-- A Real instance would make realToFrac work.
-- instance Real RealCyclotomic where
--     toRational c = toRational (toReal c)

-- | The square root of a 'Rational' number.
sqrtRat :: Rational -> RealCyclotomic
sqrtRat r
    | r >= 0  = RealCyclotomic (Cyc.sqrtRat r)
    | otherwise  = error "sqrtRational needs a nonnegative argument"

-- | Sine function with argument in degrees.
sinDeg :: Rational -> RealCyclotomic
sinDeg r = RealCyclotomic (Cyc.sinDeg r)

-- | Cosine function with argument in degrees.
cosDeg :: Rational -> RealCyclotomic
cosDeg r = RealCyclotomic (Cyc.cosDeg r)

-- | Sine function with argument in revolutions.
sinRev :: Rational -> RealCyclotomic
sinRev r = RealCyclotomic (Cyc.sinRev r)

-- | Cosine function with argument in revolutions.
cosRev :: Rational -> RealCyclotomic
cosRev r = RealCyclotomic (Cyc.cosRev r)

-- | Is the cyclotomic a rational?
isRat :: RealCyclotomic -> Bool
isRat (RealCyclotomic r) = Cyc.isRat r

-- | Return an exact rational number if possible.
toRat :: RealCyclotomic -> Maybe Rational
toRat (RealCyclotomic r) = Cyc.toRat r

-- | Export as an inexact real number.
toReal :: RealFloat a => RealCyclotomic -> a
toReal (RealCyclotomic r) = realPart (Cyc.toComplex r)

-- | The golden ratio, @(1 + √5)/2@.
goldenRatio :: RealCyclotomic
goldenRatio = (1 + sqrtRat 5) / 2

-- | Heron's formula for the area of a triangle with
--   side lengths a, b, c.
heron :: Rational        -- ^ a
      -> Rational        -- ^ b
      -> Rational        -- ^ c
      -> RealCyclotomic  -- ^ area of triangle
heron a b c
    = sqrtRat (s * (s-a) * (s-b) * (s-c))
      where
        s = (a + b + c) / 2
