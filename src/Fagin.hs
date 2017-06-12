module Fagin
(
    readModels
  , GeneModel
  , writeResult
  , writeResultAndExit
) where

import Fagin.Gff
import Fagin.Report
import Fagin.GeneModel
import System.Exit as SE
import System.IO

import qualified Data.Text as T

readModels :: T.Text -> ReportS [GeneModel]
readModels s = readGff s >>= buildModels

writeResult :: (Monoid e, ShowE e, Show o) => Report e o -> IO ()
writeResult (Pass x w n)
  =  hPutStr stdout (show x)
  >> hPutStr stderr (show3E mempty w n)
writeResult (Fail e w n)
  =  hPutStr stderr (show3E e w n)

writeResultAndExit :: (Monoid e, ShowE e, Show o) => Report e o -> IO a
writeResultAndExit (Pass x w n)
  =  writeResult (Pass x w n)
  >> SE.exitWith SE.ExitSuccess
writeResultAndExit f
  =  writeResult f
  >> SE.exitWith (SE.ExitFailure 1)
