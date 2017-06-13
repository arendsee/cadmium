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
      model_chrid  :: T.Text       -- ^ chromosome/scaffold
    , model_parent :: Maybe T.Text -- ^ parent (a 'gene') is given
    , model_id     :: T.Text       -- ^ the ID value (must be unique)
    , model_cds    :: [Interval]   -- ^ list of coding intervals
    , model_exon   :: [Interval]   -- ^ list of exon intervals
    , model_strand :: Strand       -- ^ strand [+-.?]
  }
  deriving(Show)

instance Show GeneModel where
  show gm = T.unpack model_chrid gm ++ "\t" ++ T.unpack 


type ParentId = T.Text
type EntryId = T.Text

type IdMap = MS.Map EntryId GffEntry

data Parent = Parent ParentId IntervalType deriving(Show,Eq,Ord)

requireParent :: ([Parent], GffEntry) -> ReportS ([Parent], GffEntry)
requireParent ([], g) = case gff_type g of
  Exon -> fail' $ "FeatureExpectParent: exon must have a parent\n" ++ show g
  CDS  -> fail' $ "FeatureExpectParent: CDS must have a parent\n" ++ show g
  _    -> pass' ([], g)
requireParent x = pass' x

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
      Nothing  -> fail' $ unwords
        [
            "FeatureInvalidParent: In 'Parent="
          , T.unpack p
          , "', no matchinf 'Id' found\n"
          , " - Offending line:\n"
          , show g
        ]


mapEntries :: [GffEntry] -> IdMap
mapEntries = MS.fromList . concatMap plist where
  plist g = case attrID $ gff_attr g of
    Just p  -> [(p, g)]
    Nothing -> []

toModels :: [(Parent, GffEntry)] -> ReportS [GeneModel]
toModels = sequence . map toModel . LE.groupSort where

  -- TODO check that each CDS is subsumed by an exon
  -- TODO check that no exons overlap

  toModel :: (Parent, [GffEntry]) -> ReportS GeneModel
  toModel ((Parent pid _), gs)
    = case L.group . M.catMaybes . map gff_strand $ gs of
      [(s:_)] ->  GeneModel
              <$> getChrid gs  -- model_chrid
              <*> pure Nothing -- model_parent
              <*> pure pid     -- model_id
              <*> getCDS' gs   -- model_cds
              <*> getExon' gs  -- model_exon
              <*> pure s       -- model_strand
      [] -> fail' $ "InvalidGeneModel: CDS and exon entries specify no strand\n"
      _  -> fail' $ "InvalidGeneModel: all elements of a gene model must be on the same strand"

  getChrid :: [GffEntry] -> ReportS T.Text
  getChrid gs = case L.group . map gff_seqid $ gs of
    [(s:_)] -> pass' s
    ss -> fail' $ "InvalidGeneModel: model span multiple scaffolds: " ++
                  L.intercalate ", " (map (T.unpack . head) ss)

  getCDS' :: [GffEntry] -> ReportS [Interval]
  getCDS' gs = pass' . map gff_interval . filter isCDS $ gs

  isCDS :: GffEntry -> Bool
  isCDS g = case gff_type g of
    CDS -> True
    _   -> False

  isExon :: GffEntry -> Bool
  isExon g = case gff_type g of
    Exon -> True
    _    -> False

  getExon' :: [GffEntry] -> ReportS [Interval]
  getExon' gs = pass' . map gff_interval . filter isExon $ gs


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
