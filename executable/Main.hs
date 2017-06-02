import qualified Data.Text.IO as TIO
import qualified System.IO as IO

import Fagin (readGff)

gfffile :: String
gfffile = "sample-data/test.gff"

main :: IO ()
main = do
  dat <- TIO.readFile gfffile
  case readGff dat of
    Left err -> IO.hPutStrLn IO.stderr (" *** ERROR in " ++ gfffile ++ " ***\n") >>
                IO.hPutStrLn IO.stderr (show err)
    Right gff -> IO.hPutStrLn IO.stdout (show gff)
