{-|

Perhaps the most error prone step of preparing data for Fagin is gathering
correct GFF files for each species.

Here are the lexical requirements:

 * has 9 TAB-delimited columns
 * column 1 matches the names of an entry in a genome fasta file
 * column 3 must include `mRNA`, `exon` and `CDS` labels
 * column 4 and 5 must be counting numbers where c5 > c4
 * column 7 must contain strand info
    [@+@] plus strand
    [@-@] minus stand
    [@?@] stranded but unknown
    [@.@] unstranded
 * column 9 must be a ';'-delimited string of '<tag>=<value>' entries

I also impose these semantic requirements:

 * Every exon and CDS entry must have at least one Parent, where the value
   matchs the value of an mRNA's ID
 * Every mRNA must have at least one exon
 * Every mRNA must have at least one CDS
 * Every CDS is subsumed by an exon in the same mRNA

I am systematically more permissive than required by the specification

 * column 1 can be a string of any non-TAB characters
 * column 8 (phase) can be missing for CDS (I don't use it)
 * column 6 (score) can be anything (again, I don't use it)
 * I allow a single untagged value in column 9, which I assign to ID if ID is missing

There are several pathological cases that a GFF parsers needs to be able to
deal with. These are outlined in the GFF specification:

 * single exon genes
 * polycistronic transcripts
 * genes with inteins
 * trans-spliced trancripts
 * programmed frameshift
 * operons
 * circular genomes (well, I consider this pathological)

All of these need special consideration. Currently they are not handled.

-}


module Fagin.Gff (
    readGff
  , IntervalType(..)
  , GffEntry(..)
  , Attributes(..)
  , defaultAttributes
) where

import Data.ByteString.Char8 (readInteger)

import Fagin.Prelude
import Fagin.Interval
import Fagin.Report (ReportS, pass', fail', warn')

-- | Holds the types that are currently used by Fagin. I may extend this later.
-- Since these types are required to be Sequence Ontology terms, I really ought
-- to just import the whole ontology table and allow all terms. Then, since
-- this is an ontology, I might as well port the relations between the terms.
-- This would be a time-consuming task, but may be worthwhile eventually.
data IntervalType
  = MRna
  | CDS
  | Exon
  | Gene
  | Other ByteString
  deriving(Eq,Ord,Show)

instance BShow IntervalType where
  bshow MRna      = "mRNA"
  bshow CDS       = "CDS"
  bshow Exon      = "exon"
  bshow Gene      = "gene"
  bshow (Other t) = t

-- | Attributes of a GFF entry, exactly according to the specification. The
-- data constructors nearly follow the tag names, except that they have been
-- converted to camel case, as per Haskell conventions, for example
-- "Derives_from" is converted to "DerivesFrom".
data Attributes
  = Attributes {
    -- | The unique ID for this entry. Presence in more than one GFF entry
    -- implies the entries are members of a single discontinuous feature (their
    -- types should be the same).
    attrID              :: !(Maybe ByteString)

    -- | The display name of the feature. Does not have to be unique.
    , attrName          :: !(Maybe ByteString)

    -- | List of aliases for the feature (for example locus and model ids)
    , attrAlias         :: ![ByteString]

    -- | List of parents of this feature. Indicates a part_of relationship.
    , attrParent        :: ![ByteString]

    -- | Not currently used by Fagin
    , attrTarget        :: !(Maybe ByteString)
    
    -- | Not currently used by Fagin
    , attrGap           :: !(Maybe ByteString)

    -- | Not currently used by Fagin
    , attrDerivesFrom   :: !(Maybe ByteString)

    -- | Free notes about the entry. These notes do not have to be quoted
    -- (according to the specification). Thus any special characters, which
    -- include commas, need to be encoded.
    , attrNote          :: ![ByteString]
    
    -- | A database cross reference
    , attrDbxref        :: ![ByteString]

    -- | Ontology cross reference
    , attrOntologyTerm  :: ![ByteString]

    -- | Is the sequence circular (e.g. a mitochondrial or bacterial genome)
    , attrIsCircular    :: !(Maybe Bool)

    -- | The tags defined above are all the tags with predefined meanings.
    -- Users are free to use any additional flags they desire. These tags must
    -- be lowercase, since the spec reserves uppercase tags be for future
    -- official use.
    , attrUserDefined   :: ![(ByteString, ByteString)]
  
  } deriving(Eq,Ord,Show)

defaultAttributes :: Attributes
defaultAttributes = Attributes {
      attrID            = Nothing
    , attrName          = Nothing
    , attrAlias         = []
    , attrParent        = []
    , attrTarget        = Nothing
    , attrGap           = Nothing
    , attrDerivesFrom   = Nothing
    , attrNote          = []
    , attrDbxref        = []
    , attrOntologyTerm  = []
    , attrIsCircular    = Nothing
    , attrUserDefined   = []
  }

instance BShow Attributes where
  bshow g =
    intercalate ";" $ concat
    [
        maybeAttr "ID"           (attrID           g)
      , maybeAttr "Name"         (attrName         g)
      , listAttr  "Alias"        (attrAlias        g)
      , listAttr  "Parent"       (attrParent       g)
      , maybeAttr "Target"       (attrTarget       g)
      , maybeAttr "Gap"          (attrGap          g)
      , maybeAttr "DerivesFrom"  (attrDerivesFrom  g)
      , listAttr  "Note"         (attrNote         g)
      , listAttr  "Dbxref"       (attrDbxref       g)
      , listAttr  "OntologyTerm" (attrOntologyTerm g)
      , maybe [] (\b -> if b then ["true"] else ["false"]) (attrIsCircular g)
      , map (\(k,v) -> k ++ "=" ++ v) (attrUserDefined g)
    ]
    where
    maybeAttr :: ByteString -> (Maybe ByteString) -> [ByteString]
    maybeAttr s (Just t) = [s ++ "=" ++ t]
    maybeAttr _ Nothing  = []

    listAttr :: ByteString -> [ByteString] -> [ByteString]
    listAttr _ [] = []
    listAttr s ts = [s ++ "=" ++ intercalate "," ts]


-- | Holds the data from a GFF entry that is relevant to Fagin. Some GFF
-- columns are skipped. These may be added later, but for now I don't need
-- them.
data GffEntry = GffEntry {

    -- | GFF column 1. The name of the genomic scaffold and chromosome to which
    -- the feature maps
    gff_seqid :: !ByteString

    -- | GFF column 3. The type of the feature. This must be a Sequence
    -- Ontology term or identification id.
    , gff_type :: !IntervalType

    -- | GFF column 4 and 5. The 1-based start and stop positions of this
    -- feature.
    , gff_interval :: !Interval

    {-| GFF column 7. The strand on which the interval resides. This must be
       one of the following:
        * '+' - plus sense
        * '-' - negative sense
        * '.' - strand is irrelevant
        * '?' - strand is relevant but unknown
    -}
    , gff_strand :: !(Maybe Strand)

    -- | GFF column 9. Feature attributes (see 'Attributes')
    , gff_attr :: !Attributes
  }
  deriving(Eq,Ord,Show)

-- | Construct a GffEntry from the full data of a parsed GFF line
gffEntry
  :: ByteString       -- ^ seqid
  -> ByteString       -- ^ source
  -> IntervalType -- ^ type
  -> Integer      -- ^ start
  -> Integer      -- ^ stop
  -> ByteString       -- ^ score
  -> Maybe Strand -- ^ strand
  -> ByteString       -- ^ phase
  -> Attributes   -- ^ attributes
  -> GffEntry
gffEntry seqid _ ftype start stop _ strand _ attr = 
  GffEntry {
      gff_seqid    = seqid
    , gff_type     = ftype
    , gff_interval = Interval start stop
    , gff_strand   = strand
    , gff_attr     = attr
  }

instance BShow GffEntry where
  bshow GffEntry { 
      gff_seqid    = seqid
    , gff_type     = ftype
    , gff_interval = Interval start stop
    , gff_strand   = strand
    , gff_attr     = attr
  } = intercalate "\t"
    [
        seqid
      , "."
      , bshow ftype
      , bshow start
      , bshow stop
      , "."
      , maybe "." bshow strand
      , "."
      , bshow attr
    ]


type GParser a = ByteString -> ReportS a

readInt :: GParser Integer
readInt s = case readInteger s of
  Just (x,"") -> pass' x
  _ -> fail' $ "GffParse: expected integer, found '" ++ s ++ "'"

readType :: GParser IntervalType
readType s = case s of
  ""                -> fail' "GffParse: exected <type> in column 3, found nothing"

  "gene"            -> pass' Gene
  "SO:0000704"      -> pass' Gene

  "mRNA"            -> pass' MRna
  "messenger RNA"   -> pass' MRna -- synonym
  "SO:0000234"      -> pass' MRna
  -- - this may or may not be a coding transcript
  -- - technically, mRNA is_a transcript, and a CDS or exon is only transitively
  --   a part of ta transcript.
  "transcript"      -> pass' MRna
  "SO:0000673"      -> pass' MRna

  "CDS"             -> pass' CDS
  "coding_sequence" -> pass' CDS -- synonym
  "coding sequence" -> pass' CDS -- synonym
  "SO:0000316"      -> pass' CDS

  "exon"            -> pass' Exon
  "SO:0000147"      -> pass' Exon
  -- This is slightly more specific then exon, hence the different SO id,
  -- however, it is not terribly common.
  "coding_exon"     -> pass' Exon
  "coding exon"     -> pass' Exon -- synonym
  "SO:0000195"      -> pass' Exon
  x                 -> pass' $ Other x

readStrand :: GParser (Maybe Strand)
readStrand s = case s of
  "+" -> pass' $ Just Plus
  "-" -> pass' $ Just Minus
  "." -> pass' Nothing -- The distinction between '.' and '?' is a nuance I
  "?" -> pass' Nothing -- will ignore for now
  v   -> fail' $ "GffParse: expected strand from set [+-.?], found '" ++ v

readAttributes :: GParser Attributes
readAttributes s =
  (   sequence
    . map toPair
    . map (splitSeq "=")
    . splitSeq ";"
    $ s
  ) >>= warnIfTagsRepeat >>= toAttr where

  toPair :: [ByteString] -> ReportS (ByteString, ByteString)
  toPair []             = fail' $ "GffAttribute: expected attribute (<tag>=<val>), found ''"
  toPair [tag,val]      = pass' (tag , val)
  toPair [val]          = pass' (""  , val)
  toPair fs             = fail' $ concat 
    ["GffAttribute: expected attribute (<tag>=<val>), found '", intercalate "=" fs, "'"]

  toAttr :: [(ByteString, ByteString)] -> ReportS Attributes
  toAttr attr = toAttr' attr <$> isCircular attr

  -- There is certainly a cleaner way to do this ...
  toAttr' :: [(ByteString, ByteString)] -> (Maybe Bool) -> Attributes
  toAttr' attr circular = 
    Attributes {
        attrID           = handleID attr
      , attrName         = lookup    "NAME"          attr
      , attrAlias        = maybeMany "Alias"         attr
      , attrParent       = maybeMany "Parent"        attr
      , attrTarget       = lookup    "Target"        attr
      , attrGap          = lookup    "Gap"           attr
      , attrDerivesFrom  = lookup    "Derives_from"  attr
      , attrNote         = maybeMany "Note"          attr
      , attrDbxref       = maybeMany "Dbxref"        attr
      , attrOntologyTerm = maybeMany "Ontology_term" attr
      , attrIsCircular   = circular
      , attrUserDefined  = filter isUserDefined attr
    }

  isCircular :: [(ByteString, ByteString)] -> ReportS (Maybe Bool)
  isCircular a = case lookup "Is_circular" a of
    Just "true"  -> pass' $ Just True
    Just "false" -> pass' $ Just False
    Nothing      -> pass' $ Nothing
    _            -> fail' $ "GFFAttribute: Is_circular tag must be either 'true' or 'false'"

  isUserDefined :: (ByteString, ByteString) -> Bool
  isUserDefined (t,_) = maybe False isLower (headMay t) where
    isLower = \c -> c <= 122 && c >= 97

  handleID :: [(ByteString, ByteString)] -> Maybe ByteString
  handleID a = case lookup "ID" a of
    Just x  -> Just x
    -- This interprets untagged value as an ID if no ID is provided
    -- The spec does not require this, but I add it in to handle the shit
    -- _certain_ programs throw at us.
    Nothing -> case lookup "" a of
      Just x  -> Just x
      Nothing -> Nothing

  maybeMany k attrs = case lookup k attrs of
    Just x -> splitSeq "," x
    Nothing -> []

  warnIfTagsRepeat :: [(ByteString, ByteString)] -> ReportS [(ByteString, ByteString)]
  warnIfTagsRepeat ps = case mapMaybe headMay . map (drop 1) . group . sort . map fst $ ps of
    [] -> pass' ps
    es -> pass' ps >>= warn' ("GFFAttribute: each tag may appear only once, offending tag(s): [" ++ tags ++ "]") where
      tags = intercalate ", " es


type GIFilter = (Integer,[ByteString]) -> Bool

comment' :: GIFilter
comment' (_,(x:_)) = isPrefixOf "#" x
comment' _ = True

empty' :: GIFilter
empty' (_,[]) = True
empty' _      = False

readGff :: ByteString -> ReportS [GffEntry]
readGff = 
  mapM toGff              . -- Parse GFF entry and report errors. Die on first
                            -- failed line. Most errors are highly repetitive
                            -- in GFFs, so just dying on the first failure
                            -- avoids extremely long error message. A better
                            -- solution would be to merge the similar errors.

  filter (not . empty')   . -- Filter out lines that are either empty
  filter (not . comment') . -- or start with a comment (#) character

  zip [1..]               . -- Add line numbers. This must precede filters
                            -- so line numbering in error messages is correct
                          
  map (splitSeq "\t")     . -- Break tests by line and TAB. NOTE:
  lines                     -- this allows space in fields

  where
    toGff (_, [c1,c2,c3,c4,c5,c6,c7,c8,c9])
      = gffEntry
      <$> pure           c1 -- seqid   (used as is)
      <*> pure           c2 -- source  (not used)
      <*> readType       c3 -- type    
      <*> readInt        c4 -- start   
      <*> readInt        c5 -- stop    
      <*> pure           c6 -- score   (not used)
      <*> readStrand     c7 -- strand  
      <*> pure           c8 -- phase   (not used)
      <*> readAttributes c9 -- attributes
    toGff (i,fs) = fail' $ concat ["(GFF line ", bshow i, ") Expected 9 columns, found '", bshow (length fs), "'"]
