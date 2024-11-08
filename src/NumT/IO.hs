-- | converting 'T' into Double and String
module NumT.IO where

import Control.Lens
import Control.Monad
import Data.List
import Data.Number.Flint
import Foreign.C
import GHC.IO.Handle
import NumT
import System.IO

-- | redirect stdout to replace file contents with an action's stdout,
--
-- this shouldn't be necessary but I would have to use/find the variation on
-- 'arb_printd' which writes to a buffer or file
captureOutputToFile :: FilePath -> IO () -> IO ()
captureOutputToFile filePath action = do
  withFile filePath WriteMode $ \handle -> do
    originalStdout <- hDuplicate stdout
    hDuplicateTo handle stdout
    action
    hDuplicateTo originalStdout stdout

-- | print the real part with error ie. 1.46137951490312385997144178344 +/- 5.1012e-21
printRealparts :: (?prec :: CLong) => [T] -> IO ()
printRealparts =
  mapM_
    ( \x -> do
        x <- x ?prec
        putStrLn ""
        withAcbRe x (`arb_printd` 30)
    )

-- | printRow (a,b,c) prints value a,error a,value b,error b,value c,error c
printRow :: (Each s s T T, ?prec :: CLong) => s -> IO ()
printRow row = do
  let printCols :: [IO ()]
      printCols =
        toListOf each row <&> \x -> void do
          x <- x ?prec
          withAcbRe
            x
            ( \x -> do
                arf_printd (arb_midref x) 30
                putStr ","
                putStr . show =<< mag_get_d (arb_radref x)
            )
  putStrLn ""
  sequence_ $ intersperse (putStr ",") printCols

midRad :: (?prec :: CLong) => T -> IO [Double]
midRad x = do
  x <- x ?prec
  withAcbRe_ x $ \x -> do
    CDouble mid <- arf_get_d (arb_midref x) arf_rnd_near
    CDouble rad <- mag_get_d (arb_radref x)
    return [mid, rad]
