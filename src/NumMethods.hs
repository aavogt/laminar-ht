-- | numerical methods operating on 'NumT.T'
module NumMethods where

import Control.Monad (when)
import Data.Data (Data)
import Data.Number.Flint
  ( Acb,
    acb_add_error_mag,
    acb_get_mag,
    acb_hypgeom_1f1,
    acb_overlaps,
    acb_sub,
    acb_swap,
    acb_union,
    arb_is_positive,
    mag_get_d,
    newAcb,
    withAcb,
    withMag,
    withNewAcb,
    withNewMag,
  )
import Foreign.C (CDouble)
import NumT (T, t, withAcbRe_, withAcb_, withAcb__)
import Text.Read (readMaybe)

-- | Bender Orszag 1999 8.1.16
--
-- Given xs::[T] is a list that converges to some value, @A@,
-- richardson2 n (take _N xs) evaluates the term below:
--
-- \[ \sum_{k=0}^N  \frac{A_{n+k} (n+k)^N (-1)^{k+N}}{k! (N- k)!} \]
--
-- this is a bit more general than 'richardson'. The uncertainty
-- in the result is increased by the absolute difference between
-- the result and the last element of xs.
--
-- something is wrong when I  use this, either it's the definition,
-- or it's being used wrong. The sums become very large and the error
-- bound includes 0. Maybe I the A_j are terms not sums?
richardson2 :: Integer -> [T] -> T
richardson2 n xs prec = do
  let _N = fromIntegral (length xs) - n
      fac j = fromIntegral $ product [1 .. j]
      coeffs :: [T]
      coeffs =
        [ (fromIntegral n + fromIntegral k) ^ _N * (-1) ^ (k + _N) / fac k / fac (_N - k)
          | k <- [0 .. _N]
        ]
  xEnd <- sum (zipWith (*) coeffs $ drop (fromIntegral n) xs) prec
  err <- abs (t xEnd - last xs) prec
  (errMag, _) <- withNewMag $ withAcb__ err . acb_get_mag
  withAcb_ xEnd $ withMag errMag . acb_add_error_mag

-- | Richardson extrapolation:
-- evaluate f at @h, h\/2, h\/4@
-- then infer the order and constant factor for the leading error term,
-- and increase the h\/4 term's error by this amount
richardson :: (Rational -> T) -> Rational -> T
richardson f h prec = do
  a <- f h prec
  at <- f (h / 2) prec
  att <- f (h / 4) prec
  k <- logBase 2 ((t at - t a) / (t att - t at)) prec
  -- or use att and at pairs? or att and a?
  -- why should any of them be different or preferable?
  b <- (t at - t a) / (fromRational h ** t k - fromRational (h / 2) ** t k) $ prec
  err <- t b * fromRational (h / 4) ** t k $ prec
  (errMag, _) <- withNewMag $ withAcb__ err . acb_get_mag
  withAcb_ att $ withMag errMag . acb_add_error_mag

-- | central difference for df/dx (O(h^4) error)
numdiff2 :: (T -> T) -> T -> Rational -> T
numdiff2 f x h prec = do
  x <- x prec
  x1 <- f (t x - fromRational (2 * h)) prec
  x2 <- f (t x - fromRational h) prec
  x3 <- f (t x + fromRational h) prec
  x4 <- f (t x + fromRational (2 * h)) prec
  (t x1 / 12 - (2 / 3) * t x2 + (2 / 3) * t x3 - t x4 / 12) / fromRational h $ prec

-- | central difference for df/dx (O(h^2) error)
numdiff :: (T -> T) -> T -> Rational -> T
numdiff f x h prec = do
  x <- x prec
  x1 <- f (t x - fromRational h) prec
  x2 <- f (t x + fromRational h) prec
  ((t x2 - t x1) / (2 * fromRational h)) prec

-- | wraps 'acb_hypgeom_1f1'
hypergeom1f1 :: T -> T -> T -> T
hypergeom1f1 a b z prec = do
  a <- a prec
  b <- b prec
  z <- z prec
  result <- newAcb
  let regularized = 1 -- regularized means dividing by gamma(b), but gamma(1) = 1 so it doesn't matter which this ends up being
  withAcb_ result \p -> withAcb__ a \a -> withAcb__ b \b -> withAcb__ z \z -> acb_hypgeom_1f1 p a b z regularized prec

-- | `secant f x0 x2 tol prec` finds x such that f x == 0, between x0 and x2
-- Data.Number.Flint.Arb.Calc's newton or bisect may be better
secant :: (T -> T) -> Acb -> Acb -> CDouble -> T
secant f x0 x2 tol prec = do
  y0 <- f (t x0) prec
  y2 <- f (t x2) prec
  (_, (_, (_, (_, (d, e))))) <- withNewMag \mag -> withNewAcb \w -> withAcb x0 \x0 -> withAcb x2 \x2 -> do
    acb_sub w x0 x2 prec
    acb_get_mag mag w
    d <- mag_get_d mag
    e <- acb_overlaps x0 x2
    return (d, e)
  if d < tol || e /= 0
    then fst <$> withNewAcb \r -> withAcb__ x0 \x0 -> withAcb__ x2 \x2 -> acb_union r x0 x2 prec
    else do
      x02 <- (t x0 - t x2) prec
      flip <- withAcbRe_ x02 arb_is_positive
      when (flip /= 0) $ withAcb__ x0 (withAcb__ x2 . acb_swap)
      x1 <- (t x2 - t y2 * (t x2 - t x0) / (t y2 - t y0)) prec
      y1 <- f (t x1) prec
      y0s <- withAcbRe_ y0 arb_is_positive
      y1s <- withAcbRe_ y1 arb_is_positive
      y2s <- withAcbRe_ y2 arb_is_positive
      case (y0s == y1s, y1s == y2s) of
        (True, False) -> secant f x1 x2 tol prec
        (False, True) -> secant f x0 x1 tol prec
        (_, _) -> return x1
