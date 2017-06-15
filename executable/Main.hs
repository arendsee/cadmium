import qualified System.Exit as SE
import System.IO (stderr, stdout)

import Fagin
import Fagin.Prelude

gfffile :: String
gfffile = "sample-data/z.gff3"


writeResult :: (Monoid e, ShowE e, BShow o) => Report e o -> IO ()
writeResult (Pass x w n)
  =  hPut stdout (bshow x)
  >> hPut stderr (show3E mempty w n)
writeResult (Fail e w n)
  =  hPut stderr (show3E e w n)

writeResultAndExit :: (Monoid e, ShowE e, BShow o) => Report e o -> IO a
writeResultAndExit (Pass x w n)
  =  writeResult (Pass x w n)
  >> SE.exitWith SE.ExitSuccess
writeResultAndExit f
  =  writeResult f
  >> SE.exitWith (SE.ExitFailure 1)


main :: IO ()
main = do
  text <- readFile gfffile
  writeResultAndExit
    $   readGff text
    >>= buildModels
    >>= sequence . map model2gff
