import qualified Data.Text.IO as TIO

import Fagin (readModels, writeResultAndExit)

gfffile :: String
gfffile = "sample-data/test.gff"

main :: IO a
main = do
  text <- TIO.readFile gfffile
  writeResultAndExit $ readModels text
