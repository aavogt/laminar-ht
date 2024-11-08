-- | numbers as functions of precision
module NumT where

import Control.Lens
import Control.Monad (void)
import Data.Data (Data)
import Data.Data.Lens (template)
import Data.Number.Flint
  ( Acb,
    CAcb,
    Ptr,
    acb_abs,
    acb_acos,
    acb_acosh,
    acb_add,
    acb_asin,
    acb_asinh,
    acb_atan,
    acb_atanh,
    acb_const_pi,
    acb_cos,
    acb_cosh,
    acb_div,
    acb_exp,
    acb_log,
    acb_mul,
    acb_pow,
    acb_set_arb,
    acb_set_d,
    acb_set_d_d,
    acb_sgn,
    acb_sin,
    acb_sinh,
    acb_sub,
    arb_midref,
    arb_printd,
    arb_radref,
    arf_printd,
    mag_get_d,
    newAcb,
    newArb,
    withAcb,
    withAcbRe,
    withArb,
  )
import Foreign.C (CLong)

-- $todo
--
-- consider
--
-- > Data.Number.Flint.Arb.FpWrap
-- > Data.Number.Flint.Arb.Hypgeom.arb_hypergeom_1f1 since I
-- > discard the imaginary parts without checking each time

type T = CLong -> IO Acb

t :: Acb -> T
t x _ = return x

pattern T :: (?prec :: CLong) => IO Acb -> T
pattern T x <- (\f -> f ?prec -> x)
  where
    T x _ = x

withAcb_ :: Acb -> (Ptr CAcb -> IO a) -> IO Acb
withAcb_ x f = fst <$> withAcb x f

withAcb__ :: Acb -> (Ptr CAcb -> IO ()) -> IO ()
withAcb__ x f = snd <$> withAcb x f

withAcbRe_ x f = snd <$> withAcbRe x f

-- | evaluate numbers (T) at the given precision, without changing the type
-- further evaluate will not change the result
evaluate :: (Data a, ?prec :: CLong) => a -> IO a
evaluate = template %%~ (\x -> t <$> x ?prec)

instance Data (CLong -> IO Acb) -- should not be necessary...

-- | alternative to Data.Number.Flint.Acb.ComplexField
instance Num (CLong -> IO Acb) where
  (+) x y prec = do
    x <- x prec
    y <- y prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ y \y -> withAcb__ x \x -> acb_add p x y prec

  (-) x y prec = do
    x <- x prec
    y <- y prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ y \y -> withAcb__ x \x -> acb_sub p x y prec

  (*) x y prec = do
    x <- x prec
    y <- y prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ y \y -> withAcb__ x \x -> acb_mul p x y prec

  fromInteger n prec = do
    result <- newAcb
    withAcb_ result \p -> acb_set_d_d p (fromIntegral n) 0

  abs x prec = do
    x <- x prec
    p <- newArb
    q <- newAcb
    withAcb_ q \q -> withArb p \p -> withAcb__ x \x -> do
      acb_abs p x prec
      acb_set_arb q p

  signum x prec = do
    x <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x \x -> acb_sgn p x prec

instance Fractional (CLong -> IO Acb) where
  (/) x y prec = do
    x' <- x prec
    y' <- y prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ y' \y -> withAcb__ x' \x -> acb_div p x y prec

  fromRational r prec = do
    result <- newAcb
    withAcb_ result \p -> acb_set_d p (fromRational r)

instance Floating (CLong -> IO Acb) where
  pi prec = do
    result <- newAcb
    withAcb_ result \p -> acb_const_pi p prec

  exp x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_exp p x prec

  log x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_log p x prec

  sin x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_sin p x prec

  cos x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_cos p x prec

  asin x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_asin p x prec

  acos x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_acos p x prec

  atan x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_atan p x prec

  sinh x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_sinh p x prec

  cosh x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_cosh p x prec

  asinh x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_asinh p x prec

  acosh x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_acosh p x prec

  atanh x prec = do
    x' <- x prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ x' \x -> acb_atanh p x prec

  (**) x y prec = do
    x' <- x prec
    y' <- y prec
    result <- newAcb
    withAcb_ result \p -> withAcb__ y' \y -> withAcb__ x' \x -> acb_pow p x y prec
