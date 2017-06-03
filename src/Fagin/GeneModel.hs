{-# LANGUAGE OverloadedStrings #-}

module Fagin.GeneModel (
    GeneModel(..)
  , GeneModelError
) where

import qualified Data.Map.Strict as M
import qualified Data.Text as T

import Fagin.Gff
import Fagin.Interval

data GeneModel = GeneModel {
    model_parent :: Maybe T.Text
  , model_id     :: T.Text
  , model_cds    :: [Pos]
  , model_exon   :: [Pos]
  , model_strand :: Maybe Strand
}

data GeneModelError
  = InvalidGff CGffError
  | InvalidParent
  | ExpectParent

type ParentId = T.Text
type EntryId = T.Text

type IdMap = M.Map EntryId GffEntry

type ErrIndGffEntry = Either GeneModelError (Maybe Parent, GffEntry)

data Parent = Parent ParentId IntervalType

requireParents :: ErrIndGffEntry -> ErrIndGffEntry
requireParents (Right (Nothing, g)) = case gff_type g of
  Exon -> Left ExpectParent
  CDS  -> Left ExpectParent
  _    -> Right (Nothing, g)
requireParents x = x

extractParent :: IdMap -> GffEntry -> ErrIndGffEntry
extractParent m g = case lookup "Parent" (gff_attr g) of
  Nothing -> Right $ (Nothing, g)
  Just p -> case M.lookup p m of
    Nothing -> case p of
      "-"  -> Right (Nothing, g)
      "."  -> Right (Nothing, g)
      "NA" -> Right (Nothing, g)
      _    -> Left $ InvalidParent
    Just pg -> Right $ (Just (Parent p (gff_type pg)) , g)

mapEntries :: [GffEntry] -> IdMap
mapEntries = M.fromList . concatMap plist where
  plist g = case lookup "Id" (gff_attr g) of
    Just p  -> [(p, g)]
    Nothing -> []

doerr :: [Either GeneModelError GeneModel] -> Either GeneModelError [GeneModel]
doerr = undefined

toModels :: [ErrIndGffEntry] -> Either GeneModelError GeneModel]
toModels = undefined

buildModels :: Either CGffError [GffEntry] -> Either GeneModelError [GeneModel]
buildModels (Left e) = Left (InvalidGff e)
buildModels (Right []) = Right []
buildModels (Right es) = doerr $ toModels $ requireParents $ map (extractParent $ mapEntries es) es
