module Fagin.GeneModel (
    GeneModel(..)
  , buildModels
  , model2gff
) where

import qualified Data.Map.Strict as MS
import qualified Data.List.Extra as LE

import Fagin.Prelude
import Fagin.Gff
import Fagin.Report
import Fagin.Interval

data GeneModel = GeneModel {
      model_chrid  :: !ByteString   -- ^ chromosome/scaffold
    , model_parent :: ![ByteString] -- ^ parent (a 'gene') is given
    , model_id     :: !ByteString   -- ^ the ID value (must be unique)
    , model_cds    :: ![Interval]   -- ^ list of coding intervals
    , model_exon   :: ![Interval]   -- ^ list of exon intervals
    , model_strand :: !Strand       -- ^ strand [+-.?]
  } deriving(Show)

model2gff :: GeneModel -> ReportS [GffEntry]
model2gff GeneModel {
      model_chrid  = chrid'
    , model_parent = parent'
    , model_id     = id'
    , model_cds    = cdss
    , model_exon   = (exon0:exons) -- require at least one exon
    , model_strand = strand'
  } = pass' $ [mrnaGff] ++ cdsGff ++ exonGff
  where
    mrnaGff = GffEntry {
        gff_seqid    = chrid'
      , gff_type     = MRna
      , gff_interval = foldr' (<>) exon0 (cdss ++ exons)
      , gff_strand   = Just strand'
      , gff_attr     = mrnaAttr parent'
    }

    mrnaAttr [] = [AttrID id']  
    mrnaAttr ps = [AttrID id', AttrParent ps]  

    cdsGff = map (newChild CDS) cdss
    exonGff = map (newChild Exon) (exon0:exons)

    newChild :: IntervalType -> Interval -> GffEntry
    newChild t i = GffEntry {
        gff_seqid    = chrid'
      , gff_type     = t
      , gff_interval = i
      , gff_strand   = Just strand'
      , gff_attr     = [AttrParent [id']]
    }
model2gff _ = fail' "InvalidModel: no exons"

type ParentId = ByteString
type EntryId = ByteString

type IdMap = MS.Map EntryId GffEntry

data Parent = Parent ParentId IntervalType deriving(Show,Eq,Ord)

requireParent :: ([Parent], GffEntry) -> ReportS ([Parent], GffEntry)
requireParent ([], g) = case gff_type g of
  Exon -> fail' $ "FeatureExpectParent: exon must have a parent\n" ++ bshow g
  CDS  -> fail' $ "FeatureExpectParent: CDS must have a parent\n" ++ bshow g
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

  -- [ByteString] -> [ByteString]
  -- I have to add this filter because some programs don't follow the specs.
  -- They add 'Parent' tags to things that don't have parents, like genes (i.e.
  -- "Parent=-"). They shouldn't do this, they really should know better. But
  -- I'm a lowly programer, my tools have to magically handle whatever
  -- malformed trash gets thrown at it.
  filter (/= "-") .

  -- _ -> [ByteString] -- list of Parent ids
  extractParents .

  -- GffEntry -> Attributes
  gff_attr $ g

  where
    extractParents :: [Attribute] -> [ByteString]
    extractParents attrs = case find isParent attrs of
      Just (AttrParent ps) -> ps
      _ -> []

    isParent :: Attribute -> Bool
    isParent (AttrParent _) = True
    isParent             _  = False

    getParent :: IdMap -> ByteString -> ReportS (Parent)
    getParent m' p = case MS.lookup p m' of
      Just pg  -> pass' $ Parent p (gff_type pg)
      Nothing  -> fail' $ unwords
        [
            "FeatureInvalidParent: In 'Parent="
          , p
          , "', no matchinf 'Id' found\n"
          , " - Offending line:\n"
          , bshow g
        ]


mapEntries :: [GffEntry] -> IdMap
mapEntries = MS.fromList . concatMap plist where
  plist :: GffEntry -> [(ByteString, GffEntry)]
  plist g = case find isID $ gff_attr g of
    Just (AttrID p)  -> [(p, g)]
    _ -> []

  isID :: Attribute -> Bool
  isID (AttrID _) = True
  isID _          = False

toModels :: [(Parent, GffEntry)] -> ReportS [GeneModel]
toModels = sequence . map toModel . LE.groupSort where

  -- TODO check that each CDS is subsumed by an exon
  -- TODO check that no exons overlap

  toModel :: (Parent, [GffEntry]) -> ReportS GeneModel
  toModel ((Parent pid _), gs)
    = case group . catMaybes . map gff_strand $ gs of
      [(s:_)] ->  GeneModel
              <$> getChrid gs  -- model_chrid
              <*> pure []      -- model_parent
              <*> pure pid     -- model_id
              <*> getCDS' gs   -- model_cds
              <*> getExon' gs  -- model_exon
              <*> pure s       -- model_strand
      [] -> fail' $ "InvalidGeneModel: CDS and exon entries specify no strand\n"
      _  -> fail' $ "InvalidGeneModel: all elements of a gene model must be on the same strand"

  getChrid :: [GffEntry] -> ReportS ByteString 
  getChrid gs = case group . map gff_seqid $ gs of
    [(s:_)] -> pass' s
    ss -> fail' $ "InvalidGeneModel: model span multiple scaffolds: " ++
                  unsplit ',' (mapMaybe headMay ss)

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
