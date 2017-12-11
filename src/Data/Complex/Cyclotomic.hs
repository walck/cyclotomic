{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE Trustworthy #-}

{- | 
Module      :  Data.Complex.Cyclotomic
Copyright   :  (c) Scott N. Walck 2012-2017
License     :  GPL-3 (see LICENSE)
Maintainer  :  Scott N. Walck <walck@lvc.edu>
Stability   :  experimental

The cyclotomic numbers are a subset of the complex numbers with
the following properties:
    
     1.  The cyclotomic numbers are represented exactly, enabling exact
     computations and equality comparisons.
    
     2.  The cyclotomic numbers contain the Gaussian rationals
     (complex numbers of the form 'p' + 'q' 'i' with 'p' and 'q' rational).
     As a consequence, the cyclotomic numbers are a dense subset of the
     complex numbers.
    
     3.  The cyclotomic numbers contain the square roots of all rational numbers.
    
     4.  The cyclotomic numbers form a field:  they are closed under addition, subtraction,
     multiplication, and division.
    
     5.  The cyclotomic numbers contain the sine and cosine of all rational
     multiples of pi.
    
     6.  The cyclotomic numbers can be thought of as the rational field extended
     with 'n'th roots of unity for arbitrarily large integers 'n'.

     Floating point numbers do not do well with equality comparison:

>(sqrt 2 + sqrt 3)^2 == 5 + 2 * sqrt 6
> -> False

     "Data.Complex.Cyclotomic" represents these numbers exactly, allowing equality comparison:

>(sqrtRat 2 + sqrtRat 3)^2 == 5 + 2 * sqrtRat 6
> -> True

     'Cyclotomic's can be exported as inexact complex numbers using the 'toComplex' function:

>e 6
> -> -e(3)^2
>real $ e 6
> -> 1/2
>imag $ e 6
> -> -1/2*e(12)^7 + 1/2*e(12)^11
>imag (e 6) == sqrtRat 3 / 2
> -> True
>toComplex $ e 6
> -> 0.5000000000000003 :+ 0.8660254037844384

     The algorithms for cyclotomic numbers are adapted from code by
     Martin Schoenert and Thomas Breuer in the GAP project <http://www.gap-system.org/>
     (in particular source files gap4r4\/src\/cyclotom.c and
     gap4r4\/lib\/cyclotom.gi).
-}

module Data.Complex.Cyclotomic
    ( Cyclotomic
    , i
    , e
    , sqrtInteger
    , sqrtRat
    , sinDeg
    , cosDeg
    , sinRev
    , cosRev
    , gaussianRat
    , polarRat
    , polarRatDeg
    , polarRatRev
    , conj
    , real
    , imag
    , isReal
    , isRat
    , isGaussianRat
    , toComplex
    , toReal
    , toRat
    , goldenRatio
    , dft
    , dftInv
    , rootsQuadEq
    , heron
    )
    where

import Data.List
    ( nub
    )
import Data.Ratio
    ( (%)
    , numerator
    , denominator
    )
import Data.Complex
    ( Complex(..)
    , realPart
    )
import qualified Data.Map as M
    ( Map
    , empty
    , singleton
    , lookup
    , keys
    , elems
    , size
    , fromList
    , toList
    , mapKeys
    , filter
    , insertWith
    , delete
    , map
    , unionWith
    , findWithDefault
    , fromListWith
    )
import Math.NumberTheory.Primes.Factorisation
    ( factorise
    )

-- | A cyclotomic number.
data Cyclotomic = Cyclotomic { order  :: Integer
                             , coeffs :: M.Map Integer Rational
                             } deriving (Eq)

-- | @abs@ and @signum@ are partial functions.
--   A cyclotomic number is not guaranteed to have a cyclotomic absolute value.
--   When defined, @signum c@ is the complex number with magnitude 1 that has the same argument as c;
--   @signum c = c / abs c@.
instance Num Cyclotomic where
    (+) = sumCyc
    (*) = prodCyc
    (-) c1 c2 = sumCyc c1 (aInvCyc c2)
    negate = aInvCyc
    abs = absVal
    signum = sigNum
    fromInteger 0 = zeroCyc
    fromInteger n = Cyclotomic 1 (M.singleton 0 (fromIntegral n))

instance Fractional Cyclotomic where
    recip = invCyc
    fromRational 0 = zeroCyc
    fromRational r = Cyclotomic 1 (M.singleton 0 r)

-- | The primitive @n@th root of unity.
--   For example, @'e'(4) = 'i'@ is the primitive 4th root of unity,
--   and 'e'(5) = exp(2*pi*i/5) is the primitive 5th root of unity.
--   In general, 'e' 'n' = exp(2*pi*i/'n').
e :: Integer -> Cyclotomic
e n
    | n < 1      = error "e requires a positive integer"
    | n == 1     = Cyclotomic 1 (M.singleton 0 1)
    | otherwise  = cyclotomic n $ convertToBase n (M.singleton 1 1)

instance Show Cyclotomic where
    show (Cyclotomic n mp)
        | mp == M.empty  = "0"
        | otherwise      = leadingTerm rat n ex ++ followingTerms n xs
        where ((ex,rat):xs) = M.toList mp

showBaseExp :: Integer -> Integer -> String
showBaseExp n 1  = "e(" ++ show n ++ ")"
showBaseExp n ex = "e(" ++ show n ++ ")^" ++ show ex

leadingTerm :: Rational -> Integer -> Integer -> String
leadingTerm r _ 0 = showRat r
leadingTerm r n ex
    | r == 1     = t
    | r == (-1)  = '-' : t
    | r > 0      = showRat r ++ "*" ++ t
    | r < 0      = "-" ++ showRat (abs r) ++ "*" ++ t
    | otherwise  = ""
    where t = showBaseExp n ex

followingTerms :: Integer -> [(Integer,Rational)] -> String
followingTerms _ [] = ""
followingTerms n ((ex,rat):xs) = followingTerm rat n ex ++ followingTerms n xs

followingTerm :: Rational -> Integer -> Integer -> String
followingTerm r n ex
    | r == 1     = " + " ++ t
    | r == (-1)  = " - " ++ t
    | r > 0      = " + " ++ showRat r ++ "*" ++ t
    | r < 0      = " - " ++ showRat (abs r) ++ "*" ++ t
    | otherwise  = ""
    where t = showBaseExp n ex

showRat :: Rational -> String
showRat r
    | d == 1     = show n
    | otherwise  = show n ++ "/" ++ show d
    where
      n = numerator r
      d = denominator r

-- GAP function EB from gap4r4/lib/cyclotom.gi
eb :: Integer -> Cyclotomic
eb n
    | n < 1           = error "eb needs a positive integer"
    | n `mod` 2 /= 1  = error "eb needs an odd integer"
    | n == 1          = zeroCyc
    | otherwise       = let en = e n
                        in sum [en^(k*k `mod` n) | k <- [1..(n-1) `div` 2]]

sqrt2 :: Cyclotomic
sqrt2 = e 8 - e 8 ^ (3 :: Int)

-- | The square root of an 'Integer'.
sqrtInteger :: Integer -> Cyclotomic
sqrtInteger n
    | n == 0     = zeroCyc
    | n < 0      = i * sqrtPositiveInteger (-n)
    | otherwise  = sqrtPositiveInteger n

sqrtPositiveInteger :: Integer -> Cyclotomic
sqrtPositiveInteger n
    | n < 1      = error "sqrtPositiveInteger needs a positive integer"
    | otherwise  = let factors = factorise n
                       factor = product [p^(m `div` 2) | (p,m) <- factors]
                       nn     = product [p^(m `mod` 2) | (p,m) <- factors]
                   in case nn `mod` 4 of
                        1 -> fromInteger factor * (2 * eb nn + 1)
                        2 -> fromInteger factor * sqrt2 * sqrtPositiveInteger (nn `div` 2)
                        3 -> fromInteger factor * (-i) * (2 * eb nn + 1)
                        _ -> fromInteger factor * 2 * sqrtPositiveInteger (nn `div` 4)

-- | The square root of a 'Rational' number.
sqrtRat :: Rational -> Cyclotomic
sqrtRat r = prodRatCyc (1 % fromInteger den) (sqrtInteger (numerator r * den))
    where
      den = denominator r

-- | The square root of -1.
i :: Cyclotomic
i = e 4

-- | Make a Gaussian rational; @gaussianRat p q@ is the same as @p + q * i@.
gaussianRat :: Rational -> Rational -> Cyclotomic
gaussianRat p q = fromRational p + fromRational q * i

-- | A complex number in polar form, with rational magnitude @r@ and rational angle @s@
--   of the form @r * exp(2*pi*i*s)@; @polarRat r s@ is the same as @r * e q ^ p@,
--   where @s = p/q@.  This function is the same as 'polarRatRev'.
polarRat :: Rational    -- ^ magnitude
         -> Rational    -- ^ angle, in revolutions
         -> Cyclotomic  -- ^ cyclotomic number
polarRat r s
    = let p = numerator s
          q = denominator s
      in case p >= 0 of
           True  -> fromRational r * e q ^ p
           False -> conj $ fromRational r * e q ^ (-p)

-- | A complex number in polar form, with rational magnitude and rational angle
--   in degrees.
polarRatDeg :: Rational    -- ^ magnitude
            -> Rational    -- ^ angle, in degrees
            -> Cyclotomic  -- ^ cyclotomic number
polarRatDeg r deg
    = let s = deg / 360
          p = numerator s
          q = denominator s
      in case p >= 0 of
           True  -> fromRational r * e q ^ p
           False -> conj $ fromRational r * e q ^ (-p)

-- | A complex number in polar form, with rational magnitude and rational angle
--   in revolutions.
polarRatRev :: Rational    -- ^ magnitude
            -> Rational    -- ^ angle, in revolutions
            -> Cyclotomic  -- ^ cyclotomic number
polarRatRev r s
    = let p = numerator s
          q = denominator s
      in case p >= 0 of
           True  -> fromRational r * e q ^ p
           False -> conj $ fromRational r * e q ^ (-p)

-- | Complex conjugate.
conj :: Cyclotomic -> Cyclotomic
conj (Cyclotomic n mp)
    = mkCyclotomic n (M.mapKeys (\k -> (n-k) `mod` n) mp)

-- | Real part of the cyclotomic number.
real :: Cyclotomic -> Cyclotomic
real z = (z + conj z) / 2

-- | Imaginary part of the cyclotomic number.
imag :: Cyclotomic -> Cyclotomic
imag z = (z - conj z) / (2*i)

absVal :: Cyclotomic -> Cyclotomic
absVal c = let modsq = c * conj c
           in case toRat modsq of
                Just msq -> sqrtRat msq
                Nothing  -> error "abs not available for this number"

sigNum :: Cyclotomic -> Cyclotomic
sigNum 0 = zeroCyc
sigNum c = let modsq = c * conj c
           in if isRat modsq then c / absVal c
              else error "signum not available for this number"

convertToBase :: Integer -> M.Map Integer Rational -> M.Map Integer Rational
convertToBase n mp = foldr (\(p,r) m -> replace n p r m) mp (extraneousPowers n)

removeZeros :: M.Map Integer Rational -> M.Map Integer Rational
removeZeros = M.filter (/= 0)

-- Corresponds to GAP implementation.
-- Expects that convertToBase has already been done.
cyclotomic :: Integer -> M.Map Integer Rational -> Cyclotomic
cyclotomic ord = tryReduce . tryRational . gcdReduce . Cyclotomic ord

mkCyclotomic :: Integer -> M.Map Integer Rational -> Cyclotomic
mkCyclotomic ord = cyclotomic ord . removeZeros . convertToBase ord

-- | Step 1 of cyclotomic is gcd reduction.
gcdReduce :: Cyclotomic -> Cyclotomic
gcdReduce cyc@(Cyclotomic n mp) = case gcdCyc cyc of
                                    1 -> cyc
                                    d -> Cyclotomic (n `div` d) (M.mapKeys (`div` d) mp)

gcdCyc :: Cyclotomic -> Integer
gcdCyc (Cyclotomic n mp) = gcdList (n:M.keys mp)

-- | Step 2 of cyclotomic is reduction to a rational if possible.
tryRational :: Cyclotomic -> Cyclotomic
tryRational c
    | lenCyc c == fromIntegral phi && sqfree
        = case equalCoefficients c of
            Nothing -> c
            Just r  -> fromRational $ (-1)^(nrp `mod` 2)*r
    | otherwise
        = c
    where
      (phi,nrp,sqfree) = phiNrpSqfree (order c)

-- | Compute phi(n), the number of prime factors, and test if n is square-free.
--   We do these all together for efficiency, so we only call factorise once.
phiNrpSqfree :: Integer -> (Integer,Int,Bool)
phiNrpSqfree n = (phi,nrp,sqfree)
    where
      factors = factorise n
      phi = foldr (\p n' -> n' `div` p * (p-1)) n [p | (p,_) <- factors]
      nrp = length factors
      sqfree = all (<=1) [m | (_,m) <- factors]

equalCoefficients :: Cyclotomic -> Maybe Rational
equalCoefficients (Cyclotomic _ mp)
    = case ts of
        []    -> Nothing
        (x:_) -> if equal ts then Just x else Nothing
      where
        ts = M.elems mp

lenCyc :: Cyclotomic -> Int
lenCyc (Cyclotomic _ mp) = M.size $ removeZeros mp

-- | Step 3 of cyclotomic is base reduction
tryReduce :: Cyclotomic -> Cyclotomic
tryReduce c
    = foldr reduceByPrime c squareFreeOddFactors
      where
        squareFreeOddFactors = [p | (p,m) <- factorise (order c), p > 2, m <= 1]

reduceByPrime :: Integer -> Cyclotomic -> Cyclotomic
reduceByPrime p c@(Cyclotomic n _)
    = case mapM (\r -> equalReplacements p r c) [0,p..n-p] of
        Just cfs -> Cyclotomic (n `div` p) $ removeZeros $ M.fromList $ zip [0..(n `div` p)-1] (map negate cfs)
        Nothing  -> c

equalReplacements :: Integer -> Integer -> Cyclotomic -> Maybe Rational
equalReplacements p r (Cyclotomic n mp)
    =  case [M.findWithDefault 0 k mp | k <- replacements n p r] of
         [] -> error "equalReplacements generated empty list"
         (x:xs) | equal (x:xs) -> Just x
         _ -> Nothing

replacements :: Integer -> Integer -> Integer -> [Integer]
replacements n p r = takeWhile (>= 0) [r-s,r-2*s..] ++ takeWhile (< n) [r+s,r+2*s..]
    where s = n `div` p

replace :: Integer -> Integer -> Integer -> M.Map Integer Rational -> M.Map Integer Rational
replace n p r mp = case M.lookup r mp of
                     Nothing  -> mp
                     Just rat -> foldr (\k m -> M.insertWith (+) k (-rat) m) (M.delete r mp) (replacements n p r)

includeMods :: Integer -> Integer -> Integer -> [Integer]
includeMods n q start = [start] ++ takeWhile (>= 0) [start-q,start-2*q..] ++ takeWhile (< n) [start+q,start+2*q..]

removeExps :: Integer -> Integer -> Integer -> [Integer]
removeExps n 2 q = concatMap (includeMods n q) $ map ((n `div` q) *) [q `div` 2..q-1]
removeExps n p q = concatMap (includeMods n q) $ map ((n `div` q) *) [-m..m]
    where m = (q `div` p - 1) `div` 2

pqPairs :: Integer -> [(Integer,Integer)]
pqPairs n = map (\(p,k) -> (p,p^k)) (factorise n)

extraneousPowers :: Integer -> [(Integer,Integer)]
extraneousPowers n
    | n < 1      = error "extraneousPowers needs a postive integer"
    | otherwise  = nub $ concat [[(p,r) | r <- removeExps n p q] | (p,q) <- pqPairs n]

-- | Sum of two cyclotomic numbers.
sumCyc :: Cyclotomic -> Cyclotomic -> Cyclotomic
sumCyc (Cyclotomic o1 map1) (Cyclotomic o2 map2)
    = let ord = lcm o1 o2
          m1 = ord `div` o1
          m2 = ord `div` o2
          map1' = M.mapKeys (m1*) map1
          map2' = M.mapKeys (m2*) map2
      in mkCyclotomic ord $ M.unionWith (+) map1' map2'

-- | Product of two cyclotomic numbers.
prodCyc :: Cyclotomic -> Cyclotomic -> Cyclotomic
prodCyc (Cyclotomic o1 map1) (Cyclotomic o2 map2)
    = let ord = lcm o1 o2
          m1 = ord `div` o1
          m2 = ord `div` o2
      in mkCyclotomic ord $ M.fromListWith (+)
             [((m1*e1+m2*e2) `mod` ord,c1*c2) | (e1,c1) <- M.toList map1, (e2,c2) <- M.toList map2]

-- | Product of a rational number and a cyclotomic number.
prodRatCyc :: Rational -> Cyclotomic -> Cyclotomic
prodRatCyc 0 _                   = zeroCyc
prodRatCyc r (Cyclotomic ord mp) = Cyclotomic ord $ M.map (r*) mp

-- | Additive identity.
zeroCyc :: Cyclotomic
zeroCyc = Cyclotomic 1 M.empty

-- | Additive inverse.
aInvCyc :: Cyclotomic -> Cyclotomic
aInvCyc = prodRatCyc (-1)

multiplyExponents :: Integer -> Cyclotomic -> Cyclotomic
multiplyExponents j (Cyclotomic n mp)
    | gcd j n /= 1  = error "multiplyExponents needs gcd == 1"
    | otherwise     = mkCyclotomic n (M.mapKeys (\k -> j*k `mod` n) mp)

productOfGaloisConjugates :: Cyclotomic -> Cyclotomic
productOfGaloisConjugates c
    = product [multiplyExponents j c | j <- [2..ord], gcd j ord == 1]
      where
        ord = order c

-- | Multiplicative inverse.
invCyc :: Cyclotomic -> Cyclotomic
invCyc z = case toRat (z * prod) of
             Just r  -> prodRatCyc (1 / r) prod
             Nothing -> error "invCyc:  product of Galois conjugates not rational; this is a bug, please inform package maintainer"
    where
      prod = productOfGaloisConjugates z

-- | Is the cyclotomic a real number?
isReal :: Cyclotomic -> Bool
isReal c = c == conj c

-- | Is the cyclotomic a rational?
isRat :: Cyclotomic -> Bool
isRat (Cyclotomic 1 _) = True
isRat _                = False

-- | Is the cyclotomic a Gaussian rational?
isGaussianRat :: Cyclotomic -> Bool
isGaussianRat c = isRat (real c) && isRat (imag c)

-- | Export as an inexact complex number.
toComplex :: RealFloat a => Cyclotomic -> Complex a
toComplex c = sum [fromRational r * en^p | (p,r) <- M.toList (coeffs c)]
    where en = exp (0 :+ 2*pi/n)
          n = fromIntegral (order c)

-- | Export as an inexact real number if possible.
toReal :: RealFloat a => Cyclotomic -> Maybe a
toReal c
    | isReal c   = Just $ realPart (toComplex c)
    | otherwise  = Nothing

-- | Return an exact rational number if possible.
toRat :: Cyclotomic -> Maybe Rational
toRat (Cyclotomic 1 mp)
    | mp == M.empty  = Just 0
    | otherwise      = M.lookup 0 mp
toRat _ = Nothing

-- | Sine function with argument in degrees.
sinDeg :: Rational -> Cyclotomic
sinDeg d = let n = d / 360
               nm = abs (numerator n)
               dn = denominator n
               a = e dn^nm
           in fromRational(signum d) * (a - conj a) / (2*i)

-- | Cosine function with argument in degrees.
cosDeg :: Rational -> Cyclotomic
cosDeg d = let n = d / 360
               nm = abs (numerator n)
               dn = denominator n
               a = e dn^nm
           in (a + conj a) / 2

-- | Sine function with argument in revolutions.
sinRev :: Rational -> Cyclotomic
sinRev n = let nm = abs (numerator n)
               dn = denominator n
               a = e dn^nm
           in fromRational(signum n) * (a - conj a) / (2*i)

-- | Cosine function with argument in revolutions.
cosRev :: Rational -> Cyclotomic
cosRev n = let nm = abs (numerator n)
               dn = denominator n
               a = e dn^nm
           in (a + conj a) / 2

gcdList :: [Integer] -> Integer
gcdList [] = error "gcdList called on empty list"
gcdList (n:ns) = foldr gcd n ns

equal :: Eq a => [a] -> Bool
equal [] = True
equal [_] = True
equal (x:y:ys) = x == y && equal (y:ys)

-- | The golden ratio, @(1 + √5)/2@.
goldenRatio :: Cyclotomic
goldenRatio = (1 + sqrtRat 5) / 2

-- | Discrete Fourier transform,
--   @X_k = \sum_{n=0}^{N-1} x_n \cdot e^{-i 2 \pi \frac{k}{N} n}@.
dft :: [Cyclotomic] -> [Cyclotomic]
dft cs = [sum $ zipWith (*) [conj (e m^(k*n)) | n <- [0..]] cs | k <- [0..m-1]]
          where m = fromIntegral $ length cs

-- | Inverse discrete Fourier transform,
--   @x_n = \frac{1}{N} \sum_{k=0}^{N-1} X_k \cdot e^{i 2 \pi \frac{k}{N} n}@.
dftInv :: [Cyclotomic] -> [Cyclotomic]
dftInv cs = [minv * sum (zipWith (*) [e m^(k*n) | n <- [0..]] cs) | k <- [0..m-1]]
          where m = fromIntegral $ length cs
                minv = fromRational (1 % m)

-- | Solutions to the quadratic equation
--   a x^2 + b x + c = 0.
--   Returns 'Nothing' if a == 0.
rootsQuadEq :: Rational  -- ^ a
            -> Rational  -- ^ b
            -> Rational  -- ^ c
            -> Maybe (Cyclotomic,Cyclotomic)  -- ^ roots
rootsQuadEq a b c
    | a == 0     = Nothing
    | otherwise  = Just ((-bb + sqrtDisc)/(2*aa),(-bb - sqrtDisc)/(2*aa))
    where
      aa = fromRational a
      bb = fromRational b
      sqrtDisc = sqrtRat (b*b - 4*a*c)

-- | Heron's formula for the area of a triangle with
--   side lengths a, b, c.
heron :: Rational    -- ^ a
      -> Rational    -- ^ b
      -> Rational    -- ^ c
      -> Cyclotomic  -- ^ area of triangle
heron a b c
    = sqrtRat (s * (s-a) * (s-b) * (s-c))
      where
        s = (a + b + c) / 2
