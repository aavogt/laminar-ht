{-# LANGUAGE OverloadedStrings #-}

import Args
import Control.Exception hiding (evaluate)
import Control.Lens
import Control.Monad (when, zipWithM_)
import Data.Data (Data)
import Data.Foldable
import Data.Int
import Data.List
import Database.SQLite.Simple
import Foreign.C (CLong)
import NumMethods (richardson)
import NumT
import NumT.IO
import Pipe (getTemperature, lambda, _Am)
import System.Console.CmdArgs (isNormal, whenLoud, whenNormal)

main = do
  lhts <- getArgs1
  con <- open (output_db lhts)
  createTables con
  keys <- mapM (solveToDb con) (toLaminarHT2 lhts)
  putStrLn $ intercalate "," (map show keys)

solveToDb :: Connection -> LaminarHT2 -> IO Int64
solveToDb con lht@LaminarHT2 {..} = do
  mpar <- tryInsertedId con lht
  for_ mpar \par -> do
    lambdas <- evaluate [lambda pe m (realToFrac lambdatol) | m <- [0 .. nlambda - 1]]
    whenNormal do
      putStrLn "\neigenvalues are:"
      printRealparts lambdas

    _Ams <- evaluate $ map (\x -> richardson (_Am pe x) h) lambdas
    whenNormal do
      putStrLn "\n\nAm are"
      printRealparts _Ams

    for_ (zip3 [0 :: Int ..] lambdas _Ams) \(m, l, a) -> do
      [av, apm] <- midRad a
      [lv, lpm] <- midRad l
      execute con "INSERT INTO spectrum (par, m, lambda, lambdapm, a, apm) VALUES (?, ?, ?, ?, ?, ?)" (par, m, lv, lpm, av, apm)

    temps <- getTemperature sum_n pe (zip lambdas _Ams) rzs

    for_ (zip rzs temps) \((r, z), tpm) -> do
      [t, pm] <- midRad tpm
      [r, _] <- midRad r
      [z, _] <- midRad z
      execute con "INSERT INTO temperatures (par, rho, z, temp, temp_pm) VALUES (?, ?, ?, ?, ?)" (par, r, z, t, pm)
  return (either id id mpar)
  where
    -- TODO catch errors and remove par?

    ?prec = prec

-- | @tryInsertedId con f@ runs @f con@ and returns the id of the inserted row if there was no constraint violation.
--
-- I also want to get the id of the row that already exists.
tryInsertedId :: Connection -> LaminarHT2 -> IO (Either Int64 Int64)
tryInsertedId con LaminarHT2 {..} = do
  mj <-
    (Nothing <$ insertPars con) `catch` \case
      SQLError {sqlError = ErrorConstraint} -> do
        whenLoud $ putStrLn "pars already in db, skipping..."
        lookupPars con
      e -> throwIO e
  case mj of
    Nothing -> Right <$> lastInsertRowId con
    Just j -> return (Left j)

createTables con = do
  execute_
    con
    "CREATE TABLE IF NOT EXISTS pars (id INTEGER PRIMARY KEY, pe REAL, nlambda INTEGER, lambdatol REAL, h REAL, sum_n INTEGER, prec INTEGER,\
    \ UNIQUE(pe, nlambda, lambdatol, h, sum_n, prec))"

  execute_
    con
    "CREATE TABLE IF NOT EXISTS spectrum (par INTEGER, m INTEGER, lambda REAL, lambdapm REAL, a REAL, apm REAL, \
    \ FOREIGN KEY (par) REFERENCES pars(id))"

  execute_
    con
    "CREATE TABLE IF NOT EXISTS temperatures (par INTEGER, rho REAL, z REAL, temp REAL, temp_pm REAL,\
    \ FOREIGN KEY (par) REFERENCES pars(id))"
