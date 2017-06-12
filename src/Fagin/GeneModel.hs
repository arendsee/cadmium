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

import Prelude hiding(fail)
import Fagin.Gff
import Fagin.Report
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

requireParent :: ([Parent], GffEntry) -> ReportS ([Parent], GffEntry)
requireParent ([], g) = case gff_type g of
  Exon -> fail $ ["ModelExpectParent: an exon must have a parent\n" ++ show g]
  CDS  -> fail $ ["ModelExpectParent: a CDS must have a parent\n" ++ show g]
  _    -> pass ([], g)
    
requireParent x = Right x

extractParent :: IdMap -> GffEntry -> ReportS ([Parent], GffEntry)
extractParent m g =

  -- _ -> ThrowsError ([Parent], GffEntry)
  fmap (\x -> (x, g)) .

  -- _ -> ThrowsError [Parent]
  sequence .

  -- _ -> [ThrowsError Parent]
  map (getParent m) .

  -- [Text] -> [Text]
  -- I have to add this filter because some programs don't follow the specs.
  -- They add 'Parent' tags to things that don't have parents, like genes (i.e.
  -- "Parent=-"). They shouldn't do this, they really should know better. But
  -- I'm a lowly programer, my tools have to magically handle whatever
  -- malformed trash gets thrown at it.
  filter (/= "-") .

  -- _ -> [Text] -- list of Parent ids
  attrParent .

  -- GffEntry -> Attributes
  gff_attr $ g

  where
    getParent :: IdMap -> T.Text -> ReportS (Parent)
    getParent m' p = case MS.lookup p m' of
      Just pg  -> pass' $ Parent p (gff_type pg)
      Nothing  -> fail' $ unwords [
        "ModelInvalidParent: In 'Parent=",
        T.unpack s,
        "', no matchinf 'Id' found\n",
        " - Offending line:\n",
        show g
      ]


mapEntries :: [GffEntry] -> IdMap
mapEntries = MS.fromList . concatMap plist where
  plist g = case attrID $ gff_attr g of
    Just p  -> [(p, g)]
    Nothing -> []

toModels :: [(Parent, GffEntry)] -> ReportS [GeneModel]
toModels = sequence . map toModel . LE.groupSort where

  isExon :: GffEntry -> Bool
  isExon g = case gff_type g of
    Exon -> True
    _    -> False

  isCDS :: GffEntry -> Bool
  isCDS g = case gff_type g of
    CDS -> True
    _   -> False

  -- TODO check that each CDS is subsumed by an exon
  -- TODO check that no exons overlap

  toModel :: (Parent, [GffEntry]) -> ReportS GeneModel
  toModel ((Parent pid _), gs)
    = case L.group . M.catMaybes . map gff_strand $ gs of
      [(s:_)] ->
        pass' GeneModel {
            model_parent = Nothing
          , model_id     = pid
          , model_cds    = map gff_interval . filter isCDS  $ gs
          , model_exon   = map gff_interval . filter isExon $ gs
          , model_strand = s
        }
      [] -> fail' $ "ModelStrandMissing: CDS and exon entries specify no strand\n"
      _  -> fail' $ "ModelStrandMismatch: all elements of a gene model must be on the same strand"

buildModels :: [GffEntry] -> ReportS [GeneModel]
buildModels gs
  = (sequence . map extract' $ gs) >>=
    (sequence . map requireParent) >>=
    return . concatMap mRnaParent  >>=
    toModels

  where
    extract' :: GffEntry -> ReportS ([Parent], GffEntry)
    extract' = extractParent (mapEntries gs)

    isMRna :: Parent -> Bool
    isMRna (Parent _ MRna) = True
    isMRna _               = False

    mRnaParent :: ([Parent], GffEntry) -> [(Parent, GffEntry)]
    mRnaParent (ps, g) = map (\p -> (p, g)) $ filter isMRna ps
