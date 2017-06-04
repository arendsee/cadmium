module Fagin (readModels, GeneModel) where

import Fagin.Gff
import Fagin.Error
import Fagin.GeneModel

import qualified Data.Text as T

readModels :: T.Text -> ThrowsError [GeneModel]
readModels s = readGff s >>= buildModels
