import qualified Data.Text.IO as TIO
import qualified System.IO as IO

import Fagin (readGff)

main :: IO ()
main = do
  dat <- TIO.readFile "sample-data/test.gff"
  case readGff dat of
    Left err -> IO.hPutStrLn IO.stderr "ERROR:" >>
                IO.hPutStrLn IO.stderr err
    Right gff -> IO.hPutStrLn IO.stdout (show gff)
