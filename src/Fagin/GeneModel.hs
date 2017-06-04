{-# LANGUAGE OverloadedStrings #-}

module Fagin.GeneModel (
    GeneModel(..)
  , buildModels
) where

import qualified Data.Map.Strict as M
import qualified Data.Text as T
-- import qualified Data.Either as E
-- import qualified Data.List as L
-- import qualified Data.List.Extra as LE

import Fagin.Gff
import Fagin.Error
import Fagin.Interval

data GeneModel = GeneModel {
    model_parent :: Maybe T.Text
  , model_id     :: T.Text
  , model_cds    :: [Pos]
  , model_exon   :: [Pos]
  , model_strand :: Maybe Strand
}

type ParentId = T.Text
type EntryId = T.Text

type IdMap = M.Map EntryId GffEntry

type ErrIndGffEntry = ThrowsError (Maybe Parent, GffEntry)

data Parent = Parent ParentId IntervalType

requireParents :: ErrIndGffEntry -> ErrIndGffEntry
requireParents (Right (Nothing, g)) = case gff_type g of
  Exon -> Left ModelExpectParent
  CDS  -> Left ModelExpectParent
  _    -> Right (Nothing, g)
requireParents x = x

extractType :: IdMap -> GffEntry -> ErrIndGffEntry
extractType m g = case lookup "Parent" (gff_attr g) of
  Nothing -> Right $ (Nothing, g)
  Just p -> case M.lookup p m of
    Nothing -> case p of
      "-"  -> Right (Nothing, g)
      "."  -> Right (Nothing, g)
      "NA" -> Right (Nothing, g)
      _    -> Left $ ModelInvalidParent
    Just pg -> Right $ (Just (Parent p (gff_type pg)), g)

mapEntries :: [GffEntry] -> IdMap
mapEntries = M.fromList . concatMap plist where
  plist g = case lookup "Id" (gff_attr g) of
    Just p  -> [(p, g)]
    Nothing -> []

toModels :: [ErrIndGffEntry] -> ThrowsError [GeneModel]
toModels = undefined
{- toModels es = case E.partitionEithers es of -}
{-   ([],es) ->                                -}
{-   (es,_) -> Left (foldr (<>) es)            -}


buildModels :: [GffEntry] -> ThrowsError [GeneModel]
buildModels = undefined
-- buildModels es = sequence $ toModels $ map (requireParent . (extractParent $ mapEntries es)) es
