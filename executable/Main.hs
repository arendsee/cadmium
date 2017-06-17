import Fagin
import Fagin.Prelude
import Fagin.IO

gfffile :: String
gfffile = "sample-data/micro.gff3"

main :: IO ()
main = do
  text <- readFile gfffile
  writeResultAndExit
    $   readGff text
    >>= buildModels
    >>= sequence . map model2gff
