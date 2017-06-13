import qualified Data.Text.IO as TIO

import Fagin

gfffile :: String
gfffile = "sample-data/test.gff"

main :: IO a
main = do
  text <- TIO.readFile gfffile
  -- gffs :: ReportS [[GffEntry]]
  writeResultAndExit
    $ fmap (concatMap show . concat)
    $ readModels text >>= sequence . map model2gff
   

-- model2gff :: GeneModel -> ReportS [GffEntry]
-- fmap :: (a -> b) -> m a -> m b
--
-- map mode2gff :: [GeneModel] -> [ReportS [GffEntry]]
--
-- sequence . map mode2gff :: [GeneModel] -> ReportS [[GffEntry]]
--
-- readModels :: Text -> ReportS [GeneModel]
