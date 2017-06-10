{-# LANGUAGE OverloadedStrings #-}


{-|

-- * GFF3 parsing

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
) where

import qualified Data.Text as T
import           Data.Monoid ((<>))
import qualified Data.Char as DC

import Fagin.Interval
import Fagin.Error

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
  | Other T.Text
  deriving(Show,Eq,Ord)

-- | Attributes of a GFF entry, exactly according to the specification. The
-- data constructors nearly follow the tag names, except that they have been
-- converted to camel case, as per Haskell conventions, for example
-- "Derives_from" is converted to "DerivesFrom".
data Attributes
  = Attributes {
    attrID              :: Maybe T.Text
    -- ^ The unique ID for this entry. Presence in more than one GFF entry
    -- implies the entries are members of a single discontinuous feature (their
    -- types should be the same).

    , attrName          :: Maybe T.Text
    -- ^ The display name of the feature. Does not have to be unique.

    , attrAlias         :: [T.Text]
    -- ^ List of aliases for the feature (for example locus and model ids)

    , attrParent        :: [T.Text]
    -- ^ List of parents of this feature. Indicates a part_of relationship.

    , attrTarget        :: Maybe T.Text
    -- ^ Not currently used by Fagin
    
    , attrGap           :: Maybe T.Text
    -- ^ Not currently used by Fagin

    , attrDerivesFrom   :: Maybe T.Text
    -- ^ Not currently used by Fagin

    , attrNote          :: [T.Text]
    -- ^ Free notes about the entry. These notes do not have to be quoted
    -- (according to the specification). Thus any special characters, which
    -- include commas, need to be encoded.
    
    , attrDbxref        :: [T.Text]
    -- ^ A database cross reference

    , attrOntologyTerm  :: [T.Text]
    -- ^ Ontology cross reference

    , attrIsCircular    :: Maybe Bool
    -- ^ Is the sequence circular (e.g. a mitochondrial or bacterial genome)

    , attrUserDefined   :: [(T.Text, T.Text)]

    -- ^ The tags defined above are all the tags with predefined meanings.
    -- Users are free to use any additional flags they desire. These tags must
    -- be lowercase, since the spec reserves uppercase tags be for future
    -- official use.
  
  } deriving(Show,Eq,Ord)

-- | Holds the data from a GFF entry that is relevant to Fagin. Some GFF
-- columns are skipped. These may be added later, but for now I don't need
-- them.
data GffEntry = GffEntry {
    gff_seqid    :: T.Text
    -- ^ GFF column 1. The name of the genomic scaffold and chromosome to which the feature maps

    , gff_type     :: IntervalType
    -- ^ GFF column 3. The type of the feature. This must be a Sequence Ontology term or
    -- identification id.

    , gff_interval :: Interval
    -- ^ GFF column 4 and 5. The 1-based start and stop positions of this
    -- feature.

    , gff_strand   :: Maybe Strand
    -- ^ GFF column 7. The strand on which the interval resides. This must be one of the following:
    --  * '+' - plus sense
    --  * '-' - negative sense
    --  * '.' - strand is irrelevant
    --  * '?' - strand is relevant but unknown

    , gff_attr     :: Attributes
    -- ^ GFF column 9. Feature attributes (see 'Attributes')
  }
  deriving(Show,Eq,Ord)


type GParser a = T.Text -> ThrowsError a

readInt :: GParser Integer
readInt s = case reads (T.unpack s) :: [(Integer,String)] of
  [(x,"")] -> Right x
  _        -> Left $ GffExpectInteger s

readType :: GParser IntervalType
readType s = case s of
  ""                -> Left  GffNoType

  "gene"            -> Right Gene
  "SO:0000704"      -> Right Gene

  "mRNA"            -> Right MRna
  "messenger RNA"   -> Right MRna -- synonym
  "SO:0000234"      -> Right MRna
  -- * this may or may not be a coding transcript
  -- * technically, mRNA is_a transcript, and a CDS or exon is only transitively
  --   a part of ta transcript.
  "transcript"      -> Right MRna
  "SO:0000673"      -> Right MRna

  "CDS"             -> Right CDS
  "coding_sequence" -> Right CDS -- synonym
  "coding sequence" -> Right CDS -- synonym
  "SO:0000316"      -> Right CDS

  "exon"            -> Right Exon
  "SO:0000147"      -> Right Exon
  -- This is slightly more specific then exon, hence the different SO id,
  -- however, it is not terribly common.
  "coding_exon"     -> Right Exon
  "coding exon"     -> Right Exon -- synonym
  "SO:0000195"      -> Right Exon
  x                 -> Right $ Other x

readStrand :: GParser (Maybe Strand)
readStrand s = case s of
  "+" -> Right $ Just Plus
  "-" -> Right $ Just Minus
  "." -> Right Nothing
  v   -> Left $ GffExpectStrand v

readAttributes :: GParser Attributes
readAttributes s = (sequence . map toPair . map (T.splitOn "=") . T.splitOn ";" $ s) >>= toAttr where

  toPair :: [T.Text] -> ThrowsError (T.Text, T.Text)
  toPair []             = Left $ GffExpectAttribute ""
  toPair [tag,val]      = Right (tag , val)
  toPair [val]          = Right (""  , val)
  toPair fs             = Left $ GffExpectAttribute $ T.intercalate "=" fs

  toAttr :: [(T.Text, T.Text)] -> ThrowsError Attributes
  toAttr attr = case isCircular attr of
    Right circular ->
      Right $ Attributes {
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
    Left msg -> Left msg
    where
      isCircular :: [(T.Text, T.Text)] -> ThrowsError (Maybe Bool)
      isCircular a = case lookup "Is_circular" a of
        Just "true"  -> Right $ Just True
        Just "false" -> Right $ Just False
        Nothing      -> Right $ Nothing
        _            -> Left $ GFFAttributeError "Is_circular tag must be either 'true' or 'false'"

  isUserDefined :: (T.Text, T.Text) -> Bool
  isUserDefined (t,_) = T.length t > 0 && DC.isLower (T.head t)

  handleID :: [(T.Text, T.Text)] -> Maybe T.Text
  handleID a = case lookup "ID" a of
    Just x  -> Just x
    -- This interprets untagged value as an ID if no ID is provided
    -- The spec does not require this, but I add it in to handle the shit
    -- *certain* programs throw at us.
    Nothing -> case lookup "" a of
      Just x  -> Just x
      Nothing -> Nothing

  maybeMany k attrs = case lookup k attrs of
    Just x -> T.splitOn "," x
    Nothing -> []


type GIFilter = (Integer,[T.Text]) -> Bool

comment :: GIFilter
comment (_,(x:_)) = T.isPrefixOf (T.singleton '#') x
comment _ = True

empty :: GIFilter
empty (_,[]) = True
empty _      = False

readGff :: T.Text -> ThrowsError [GffEntry]
readGff =
  sequence               . -- Merge errors, die on first failure

  toGff                  . -- Parse GFF entry and report errors

  filter (not . empty)   . -- Filter out lines that are either empty
  filter (not . comment) . -- or start with a comment (#) character

  zip [1..]              . -- Add line numbers. This must precede filters
                           -- so line numbering in error messages is correct

  map (T.splitOn "\t")   . -- Break tests by line and TAB. NOTE:
  T.lines                  -- this allows space in fields
  where
    toGff :: [(Integer,[T.Text])] -> [ThrowsError GffEntry]
    toGff ((i, [chr, _, typ, a, b, _, str, _, attr]):xs)
      = case (readType typ, readInt a, readInt b, readStrand str, readAttributes attr) of
        -- If everything is good, make an entry
        (Right typ', Right a', Right b', Right str', Right attr') ->
          [ Right $
              GffEntry {
                  gff_seqid    = chr
                , gff_type     = typ'
                , gff_interval = Interval a' b'
                , gff_strand   = str'
                , gff_attr     = attr'
              }
          ] ++ toGff xs
        -- else, join the errors and annotate them with the line number
        (typ', a', b', str', attr') -> [ Left (GffLineError i err) ] where
          l :: ThrowsError a -> FaginError
          l (Right _) = mempty
          l (Left  e) = e
          err = l typ' <> l a' <> l b' <> l str' <> l attr'
    toGff [] = []
    toGff ((i,fs):_) = [ Left $ GffLineError i (GffInvalidRowNumber fs) ]
