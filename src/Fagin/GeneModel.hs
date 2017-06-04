{-# LANGUAGE OverloadedStrings #-}

module Fagin.GeneModel (
    GeneModel(..)
  , buildModels
) where

import qualified Data.Map.Strict as MS
import qualified Data.Text as T
import qualified Data.List as L
import qualified Data.List.Extra as LE
import qualified Data.Maybe as M

import Fagin.Gff
import Fagin.Error
import Fagin.Interval

data GeneModel = GeneModel {
      model_parent :: Maybe T.Text
    , model_id     :: T.Text
    , model_cds    :: [Interval]
    , model_exon   :: [Interval]
    , model_strand :: Strand
  }
  deriving(Show)

type ParentId = T.Text
type EntryId = T.Text

type IdMap = MS.Map EntryId GffEntry

data Parent = Parent ParentId IntervalType deriving(Show,Eq,Ord)

requireParent :: (Maybe Parent, GffEntry) -> ThrowsError (Maybe Parent, GffEntry)
requireParent (Nothing, g) = case gff_type g of
  Exon -> Left ModelExpectParent
  CDS  -> Left ModelExpectParent
  _    -> Right (Nothing, g)
requireParent x = Right x

extractParent :: IdMap -> GffEntry -> ThrowsError (Maybe Parent, GffEntry)
extractParent m g = case lookup "Parent" (gff_attr g) of
  Nothing -> Right $ (Nothing, g)
  Just p -> case MS.lookup p m of
    Nothing -> case p of
      "-"  -> Right (Nothing, g)
      "."  -> Right (Nothing, g)
      "NA" -> Right (Nothing, g)
      _    -> Left $ ModelInvalidParent
    Just pg -> Right $ (Just (Parent p (gff_type pg)), g)

mapEntries :: [GffEntry] -> IdMap
mapEntries = MS.fromList . concatMap plist where
  plist g = case lookup "Id" (gff_attr g) of
    Just p  -> [(p, g)]
    Nothing -> []

toModels :: [(Parent, GffEntry)] -> ThrowsError [GeneModel]
toModels = sequence . map toModel . LE.groupSort where

  isExon :: GffEntry -> Bool
  isExon g = case gff_type g of
    Exon -> True
    _    -> False

  isCDS :: GffEntry -> Bool
  isCDS g = case gff_type g of
    CDS -> True
    _   -> False

  toModel :: (Parent, [GffEntry]) -> ThrowsError GeneModel
  toModel ((Parent pid _), gs)
    = case L.group . M.catMaybes . map gff_strand $ gs of
      [(s:_)] ->
        Right $ GeneModel {
            model_parent = Nothing
          , model_id     = pid
          , model_cds    = map gff_interval . filter isCDS  $ gs
          , model_exon   = map gff_interval . filter isExon $ gs
          , model_strand = s
        }
      [] -> Left $ ModelStrandMissing
      _  -> Left $ ModelStrandMismatch

buildModels :: [GffEntry] -> ThrowsError [GeneModel]
buildModels gs
  = (sequence . map extract' $ gs) >>=
    (sequence . map requireParent) >>=
    return . removeOrphans         >>=
    toModels

  where
    extract' :: GffEntry -> ThrowsError (Maybe Parent, GffEntry)
    extract' = extractParent (mapEntries gs)

    -- remove entries without a defined parent
    removeOrphans :: [(Maybe Parent, GffEntry)] -> [(Parent, GffEntry)]
    removeOrphans ((Nothing, _):xs) = removeOrphans xs
    removeOrphans ((Just p, g):xs)  = [(p,g)] ++ removeOrphans xs
    removeOrphans [] = []
