{-# HLINT ignore "Parenthesize unary negation" #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TupleSections #-}
{-# OPTIONS_GHC -Wno-unrecognised-pragmas #-}

-- | Following Polyanin2002 section 3.5 or Lehktmakher1971
--
-- Given steady laminar flow in a cylindrical pipe, with a
-- stepwise increase in wall temperature at z=0,
-- these functions calculate the temperature field at any
-- point in the fluid after the increase (z>0).
--
--
-- \[
--  \textbf{Pe} \left(1 - \rho^2\right) T_{,z} = T_{,\rho\rho} + \rho^{-1} T_{,\rho} + T_{,zz}
-- \]
-- where dimensionless variables above are defined in terms of dimensional
-- quantities pipe radius \(R\), temperature \(T_*\), upstream wall temperature
-- \(T_1\), downstream wall temperature \(T_2\), thermal diffusivity
-- \(\alpha\), and maximum centerline velocity \(U_{\text{max}}\) as follows:
--
-- \[ \rho = \frac{r}{R} \] \[ z = \frac{Z}{R} \] \[ T =
-- \frac{T_*-T_1}{T_2-T_1} \] \[ \textbf{Pe} = R U_{\text{max}} / \alpha \]
module Pipe where

import Control.Lens
import Data.Data (Data)
import Data.Function (on)
import Data.IntMap (IntMap)
import qualified Data.IntMap as M
import Data.List
import Data.Traversable (for)
import qualified Data.Vector.Storable as V
import Foreign.C (CInt, CLong)
import Foreign.C.Types (CDouble)
import GHC.Exts
import Linear (V3 (V3))
import NumMethods
import NumT (T, t)
import Numeric.LinearAlgebra (cgSolve)
import Numeric.LinearAlgebra.Devel

-- | finite differences, rhos and zs can be unevenly spaced, but the grid is rectangular
-- done but not tested/compared with the analytical solution. This one should be better near
-- the step change in wall temperature at least.
getTemperatureFD :: Double -> V.Vector Double -> V.Vector Double -> V.Vector Double -> V.Vector (V3 Double)
getTemperatureFD pe rhos zs tWall | V.length zs == V.length tWall = V.imap addRZ x
  where
    addRZ m x =
      let (i, j) = m ^. Pipe.mij nr
       in V3 (rhos V.! i) (zs V.! j) x
    x = cgSolve False m (V.fromList b)
    m = fromCSR $ mkCSR $ concat ms

    nr = V.length rhos
    nz = V.length zs

    (ms, b) =
      unzip $
        zipWith
          (\(FDV b m) k -> (M.assocs m <&> _1 %~ (k,), b))
          (centerline ++ wall ++ zEnds ++ z0s ++ insideK)
          [0 ..]

    mij = Pipe.mij nr

    -- 4 boundary edges
    -- 2 at z=0 and z=end
    zEnds =
      [ let m = mij # (i, nz - 1)
         in fdv (V.last tWall) [(m, 1)]
        | (i, rho) <- zip [0 ..] (toList rhos)
      ]

    z0s =
      [ let m = mij # (i, 0)
         in fdv (V.head tWall) [(m, 1)]
        | (i, rho) <- zip [0 ..] (toList rhos)
      ]

    -- at rho=0 and rho=end
    wall =
      [ let m = mij # (nr - 1, j)
         in fdv wallT [(m, 1)]
        | (j, z, wallT) <- zip3 [0 ..] (toList zs) (toList tWall)
      ]

    centerline =
      [ let m01 = mij # (0, j)
            m11 = mij # (1, j)
         in fdv0 [(m01, 1), (m11, -1)]
        | (j, z) <- zip [0 ..] $ toList zs
      ]

    -- stencil points inside
    insideK =
      [ -- each equation k gets it's own row
        let -- each grid point in the temperature field gets it's own column,
            -- express the (i,j) index in terms of a single number m,
            -- where m11 is the center of the stencil
            m01 = mij # (i - 1, j)
            m11 = mij # (i, j)
            m21 = mij # (i + 1, j)
            m12 = mij # (i, j - 1)
            m10 = mij # (i, j + 1)

            -- the step sizes in rho and z
            h01 = rhos V.! i - rhos V.! (i - 1)
            h21 = rhos V.! (i + 1) - rhos V.! i
            h10 = zs V.! j - zs V.! (j - 1)
            h12 = zs V.! (j + 1) - zs V.! j

            hr = h01 + h21
            hz = h10 + h12

            -- first derivative in rho
            d10 = fdv0 [(m01, -1 / hr), (m21, 1 / hr)]

            -- first derivative in z
            d01 = fdv0 [(m10, -1 / hz), (m12, 1 / hz)]

            -- denominator in the second derivative in rho
            -- see stencil.wxmx for the derivation
            hrr = (h01 * h21 ^ 2 + h01 ^ 2 * h21) / 2

            d20 = fdv0 [(m01, h21 / hrr), (m11, -hr / hrr), (m21, h01 / hrr)]
            -- second derivative in z
            hzz = (h10 * h12 ^ 2 + h10 ^ 2 * h12) / 2
            d02 = fdv0 [(m01, h12 / hzz), (m11, -hz / hzz), (m12, h10 / hzz)]

            negPeRho = -pe * (1 - rho ^ 2)
         in fmap (negPeRho *) d01 + d20 + fmap (/ rho) d10 + d02
        | (i, rho) <- zip [1 ..] $ init $ drop 1 $ toList rhos,
          (j, z) <- zip [1 ..] $ init $ drop 1 $ toList zs
      ]

-- | terms in a finite difference equation
--
-- [(Int, a)] or IntMap a would suffice. I only have a single k at a time.
--
-- should it also contribute to the right hand side?
-- or the RHS can be part of the matrix?
data FDV a = FDV
  { fdvRhs :: a,
    unFDV :: IntMap a
  }
  deriving (Show, Functor)

-- | with RHS=0
fdv0 xs = FDV 0 $ M.fromListWith (+) xs

-- | specified right hand side
fdv n xs = FDV n $ M.fromListWith (+) xs

instance (Num a) => Num (FDV a) where
  FDV a xs + FDV b ys = FDV (a + b) $ M.unionWith (+) xs ys
  FDV a xs - FDV b ys = FDV (a - b) $ M.unionWith (-) xs ys
  FDV a xs * FDV b ys = FDV (a * b) $ M.unionWith (*) xs ys

mij :: Int -> Iso' Int (Int, Int)
mij w = iso (`divMod` w) (\(i, j) -> i * w + j)

-- should I get the hz hrho step sizes more efficiently?
-- That is, zipWith (-) zs (tail zs)
-- But the ends don't have a next element.
-- There is a leftward and a rightward step size,
-- Also diagonal?
-- or just look them up as needed?

-- | \( T = 1 - \sum_{m=0}^{\infty} A_m \exp\left(-\frac{\lambda_m^2 z}{\textbf{Pe}} - \frac{\lambda_m \rho^2}{2}\right) {}_1F_1\left(0.5 - \frac{\lambda_m}{4} - \frac{\lambda ^3}{4 \textbf{Pe}^2}, 1, \lambda_m \rho^2\right) \)
-- (eq. 3.5.8 and 3.5.14)
getTemperature ::
  Integer ->
  -- | \(\textbf{Pe}\)
  T ->
  -- | \(\lambda_m, A_m\)
  [(T, T)] ->
  -- | \(\rho_i, z_i\)
  [(T, T)] ->
  -- | \(T_i\)
  IO [T]
getTemperature n pe lambdaAms rhoZs = do
  fs :: [(T, T) -> T] <- for lambdaAms \(x, _Am) -> do
    let am = 0.5 - x / 4 - x ^ 3 / 4 / (pe ^ 2)
    return \(rho, z) prec -> do
      let zterm = -(x ^ 2) * z / pe
          rhoterm = -x * rho ^ 2 / 2
      _Am * exp (zterm + rhoterm) * hypergeom1f1 am 1 (x * rho ^ 2) $ prec
  return [1 - sum (map ($ rhoZ) fs) | rhoZ <- rhoZs]

-- | This starts with an approximation \( \lambda_m \approx 4m  + 2.7 \) (eq. 3.5.19),
-- and improves it with 'secant' to satisfy
--
-- \( 0 = {}_1F_1\left(0.5 - \frac{\lambda_m}{4} - \frac{\lambda_m ^3}{4 \textbf{Pe}^2}, 1, \lambda_m \right) \) (eq 3.5.15)
lambda ::
  -- | \(\textbf{Pe}\)
  T ->
  -- | which eigenvalue m
  Int ->
  -- | 'secant' xtol means the result uncertainty will not be less than \(\pm \texttt{xtol}/2\).
  CDouble ->
  -- | eigenvalue \(\lambda_m\)
  T
lambda pe m xtol prec = do
  pe <- pe prec
  x1 <- ((fromIntegral m * 4 + 2.7) * 1.05) prec
  x2 <- ((fromIntegral m * 4 + 2.7) * 0.95) prec
  secant
    ( \x prec -> do
        x <- x prec
        let am = 0.5 - t x / 4 - t x ^ 3 / 4 / t pe
        hypergeom1f1 am 1 (t x) prec
    )
    x1
    x2
    xtol
    prec

-- | Equation 3.5.13
--
-- \[ A_m = \frac{-2}{\lambda_m \frac{d}{d\lambda} \left[ \exp\left(-\lambda/2\right) {}_1F_1\left(0.5 - \frac{\lambda}{4} - \frac{\lambda ^3}{4 \textbf{Pe}^2}, 1, \lambda\right) \right] \Bigg|_{\lambda = \lambda_m}} \]
_Am ::
  -- | \(\textbf{Pe}\)
  T ->
  -- | \(\lambda_m\)
  T ->
  -- | @h@ step size for 'numdiff2'
  Rational ->
  T
_Am pe lambda h = (-2) / lambda / numdiff2 (amDen pe) lambda h

-- | 3.5.14 denominator for '_Am'
amDen :: T -> T -> T
amDen pe lambda prec = do
  lambda <- lambda prec
  am <- (0.5 - t lambda / 4 - t lambda ^ 3 / 4 / pe) prec
  exp (-(t lambda / 2)) * hypergeom1f1 (t am) 1 (t lambda) $ prec
