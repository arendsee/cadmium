import Criterion.Main
import qualified Control.Monad as CM
import qualified Prelude as P
import qualified Data.List as DL

import Fagin.Prelude
import Fagin.Gff
import Fagin.GeneModel

data Columns = Columns {
      coltype   :: [ByteString]
    , colstart  :: [ByteString]
    , colstop   :: [ByteString]
    , colstrand :: [ByteString]
    , colattr   :: [ByteString]
  }

getColumns :: ByteString -> Columns
getColumns rawGff
  = case   DL.transpose
         . map (split '\t')
         . filter (not . isPrefixOf "#")
         . split '\n'
         $ rawGff
    of
    [_, _, types, starts, stops, _, strands, _, attrs] -> Columns types starts stops strands attrs
    _ -> P.error "Aww fuck, that sample GFF is shit"

main :: IO ()
main = do
  rawGff <- readFile "sample-data/short.gff3"
  let objGff = readGff rawGff
  let objMod = objGff >>= buildModels
  let outgff = objMod >>= sequence . map model2gff
  let cols   = getColumns rawGff

  defaultMain [
        bgroup "fieldParsers" [
          bench "map readType"       $ nf (map readType''       ) (coltype   cols)
        , bench "map readInt"        $ nf (map readInt''        ) (colstart  cols)
        , bench "map readInt"        $ nf (map readInt''        ) (colstop   cols)
        , bench "map readStrand"     $ nf (map readStrand''     ) (colstrand cols)
        , bench "map readAttributes" $ nf (map readAttributes'' ) (colattr   cols)
      ]
      , bgroup "500kb-sample" [
          bench "readFile" $ nfIO (readFile "sample-data/short.gff3")
        , bench "readGff" $ nf readGff rawGff
        , bench "models >>= sequence . map model2gff" $ nf (CM.liftM $ sequence . map model2gff) objMod
        , bench "outgff >>= map bshow" $ nf (CM.liftM $ map bshow) outgff
      ]
    ]
