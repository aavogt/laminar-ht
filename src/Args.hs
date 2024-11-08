{-# LANGUAGE OverloadedStrings #-}

module Args where

import Control.Lens
import Data.Data (Data)
import Data.Int
import Data.Maybe (fromMaybe)
import Database.SQLite.Simple (Connection, Only (..), SQLData, execute, execute_, query)
import Database.SQLite.Simple.ToField (ToField (toField))
import Foreign.C (CLong)
import NumT
import NumT.IO
import System.Console.CmdArgs
import Text.Read (readMaybe)
import Text.Regex.Applicative
import Text.Regex.Applicative.Common

getArgs1 = fmap parse <$> cmdArgs laminarHT

data LaminarHT a = LaminarHT
  { pes, rhos, zs, nlambdas, lambdatols, hs, sum_ns, acb_precs :: a,
    output_db :: String
  }
  deriving (Data, Functor, Foldable, Traversable)

laminarHT =
  LaminarHT
    { pes = "200,500" &= help "Peclet number: 200,500",
      rhos = "0,0.05..1" &= help "radius rho: 0,0.05..1",
      zs = "0,2..100" &= help "distance z: 0,2..100",
      nlambdas = "6" &= help "number of eigenvalues: 6",
      lambdatols = "1e-20" &= help "xtolerance for eigenvalue secant method: 1e-20",
      hs = "1e-2" &= help "finite difference step size h: 1e-4",
      sum_ns = "3" &= help "number of terms in Richardson extrapolation of the eigenvalue sum (currently unused): 3",
      acb_precs = "2048" &= help "Flint2 acb precision: 2048",
      output_db = "laminarHT.db" &= help "sqlite output file: laminarHT.db"
    }
    &= verbosity
    &= help
      "Solve the laminar flow heat transfer problem with constant wall temperatures,\
      \ and store the results in an sqlite database.\
      \ Flags can have multiple values with first,step..last or a,b,c,d interpreted as usual,\
      \ and all possible combinations (cartesian product) are solved.\
      \ Results already in the database are not recomputed."
    &= details
      [ "The schema is",
        "    CREATE TABLE pars (id INTEGER PRIMARY KEY, pe REAL, nlambda INTEGER, lambdatol REAL, h REAL, sum_n INTEGER, prec INTEGER, UNIQUE(pe, nlambda, lambdatol, h, sum_n, prec));",
        "",
        "    CREATE TABLE spectrum (par INTEGER, m INTEGER, lambda REAL, lambdapm REAL, a REAL, apm REAL,  FOREIGN KEY (par) REFERENCES pars(id));",
        "",
        "    CREATE TABLE temperatures (par INTEGER, rho REAL, z REAL, temp REAL, temp_pm REAL, FOREIGN KEY (par) REFERENCES pars(id)); ",
        "",
        " `pars` stores command line arguments other than rhos and zs",
        " `spectrum` stores intermediate values",
        " `temperatures` which stores the T(rho, z)",
        "",
        " A trailing \"pm\" means it is the radius of the interval"
      ]

-- | with regex-applicative, turn a string with commas and .. into calls to enumFrom, enumFromThenTo etc.
parse :: String -> [Rational]
parse =
  maybe [] fst . findLongestPrefix do
    start <- float
    next <- many $ sym ',' *> float
    upto <- optional $ string ".." *> float
    pure $ case (start, next, upto) of
      (a, [b], Just c) -> enumFromThenTo a b c
      (a, [], Just b) -> enumFromTo a b
      (a, bs, Nothing) -> a : bs
  where
    float :: RE Char Rational
    float = do
      a <- signed decimal
      bblen <- optional $ sym '.' *> withMatched decimal
      e <- optional $ psym (`elem` ("Ee" :: String)) *> signed decimal
      pure $
        let dec = maybe 0 (\(b, blen) -> b / 10 ^ length blen) bblen
         in (a + dec) * maybe 1 (10 ^^) e

data LaminarHT2 = LaminarHT2
  { pe :: T,
    rzs :: [(T, T)],
    nlambda :: Int,
    lambdatol :: Rational,
    h :: Rational,
    sum_n :: Integer,
    prec :: CLong,
    insertPars :: Connection -> IO (),
    lookupPars :: Connection -> IO (Maybe Int64)
  }

toLaminarHT2 :: LaminarHT [Rational] -> [LaminarHT2]
toLaminarHT2 LaminarHT {..} =
  [ LaminarHT2
      { pe = fromRational pe,
        rzs = liftA2 (\a b -> (fromRational a, fromRational b)) rhos zs,
        nlambda = round nlambda,
        lambdatol = lambdatol,
        h = h,
        sum_n = round sum_n,
        prec = round prec,
        lookupPars = \con ->
          query
            con
            "SELECT id FROM pars WHERE pe = ? AND nlambda = ? AND lambdatol = ? AND h = ? AND sum_n = ? AND prec = ?"
            [ toField @Int (round pe),
              toField @Int (round nlambda),
              toField @Double (realToFrac lambdatol),
              toField @Double (realToFrac h),
              toField @Int (round sum_n),
              toField @Int (round prec)
            ]
            <&> \case
              [Only i] -> Just i
              [] -> Nothing,
        insertPars = \con ->
          execute
            con
            "INSERT INTO pars (pe, nlambda, lambdatol, h, sum_n, prec) VALUES (?, ?, ?, ?, ?, ?)"
            [ toField @Int (round pe),
              toField @Int (round nlambda),
              toField @Double (realToFrac lambdatol),
              toField @Double (realToFrac h),
              toField @Int (round sum_n),
              toField @Int (round prec)
            ]
      }
    | pe <- pes,
      nlambda <- nlambdas,
      lambdatol <- lambdatols,
      h <- hs,
      sum_n <- sum_ns,
      prec <- acb_precs
  ]
