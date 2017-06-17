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
  , Attribute(..)
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
  | Other ShortByteString
  deriving(Eq,Ord,Show)

instance BShow IntervalType where
  bshow MRna      = "mRNA"
  bshow CDS       = "CDS"
  bshow Exon      = "exon"
  bshow Gene      = "gene"
  bshow (Other t) = fromShort t

-- | Attributes of a GFF entry, exactly according to the specification. The
-- data constructors nearly follow the tag names, except that they have been
-- converted to camel case, as per Haskell conventions, for example
-- "Derives_from" is converted to "DerivesFrom".
data Attribute
    -- | The unique ID for this entry. Presence in more than one GFF entry
    -- implies the entries are members of a single discontinuous feature (their
    -- types should be the same).
    = AttrID ShortByteString

    -- | The display name of the feature. Does not have to be unique.
    | AttrName ShortByteString

    -- | List of aliases for the feature (for example locus and model ids)
    | AttrAlias [ShortByteString]

    -- | List of parents of this feature. Indicates a part_of relationship.
    | AttrParent [ShortByteString]

    -- | Not currently used by Fagin
    | AttrTarget ShortByteString
    
    -- | Not currently used by Fagin
    | AttrGap ShortByteString

    -- | Not currently used by Fagin
    | AttrDerivesFrom ShortByteString

    -- | Free notes about the entry. These notes do not have to be quoted
    -- (according to the specification). Thus any special characters, which
    -- include commas, need to be encoded.
    | AttrNote [ShortByteString]

    -- | A database cross reference
    | AttrDbxref [ShortByteString]

    -- | Ontology cross reference
    | AttrOntologyTerm [ShortByteString]

    -- | Is the sequence circular (e.g. a mitochondrial or bacterial genome)
    | AttrIsCircular Bool

    -- | The tags defined above are all the tags with predefined meanings.
    -- Users are free to use any additional flags they desire. These tags must
    -- be lowercase, since the spec reserves uppercase tags be for future
    -- official use.
    | AttrUserDefined ShortByteString ShortByteString

    -- | According to the spec, all entries in the attribute column should be
    -- wrapped in tag-value pairs. However, it is common for programs to add
    -- untagged values, which sometimes function as ID.
    | AttrUntagged ShortByteString
    deriving(Eq,Ord,Show)

instance BShow Attribute where
  bshow (AttrID           s     ) = "ID="            ++ fromShort s
  bshow (AttrName         s     ) = "Name="          ++ fromShort s
  bshow (AttrAlias        ss    ) = "Alias="         ++ unsplit ',' (map fromShort ss)
  bshow (AttrParent       ss    ) = "Parent="        ++ unsplit ',' (map fromShort ss)
  bshow (AttrTarget       s     ) = "Target="        ++ fromShort s
  bshow (AttrGap          s     ) = "Gap="           ++ fromShort s
  bshow (AttrDerivesFrom  s     ) = "Derives_from="  ++ fromShort s
  bshow (AttrNote         ss    ) = "Note="          ++ unsplit ',' (map fromShort ss)
  bshow (AttrDbxref       ss    ) = "Dbxref="        ++ unsplit ',' (map fromShort ss)
  bshow (AttrOntologyTerm ss    ) = "Ontology_term=" ++ unsplit ',' (map fromShort ss)
  bshow (AttrIsCircular   True  ) = "Is_circular=true"
  bshow (AttrIsCircular   False ) = "Is_circular=false"
  bshow (AttrUserDefined  t v   ) = fromShort t ++ "=" ++ fromShort v
  bshow (AttrUntagged     s     ) = "Untagged=" ++ fromShort s


-- | Holds the data from a GFF entry that is relevant to Fagin. Some GFF
-- columns are skipped. These may be added later, but for now I don't need
-- them.
data GffEntry = GffEntry {

    -- | GFF column 1. The name of the genomic scaffold and chromosome to which
    -- the feature maps
    gff_seqid :: !ShortByteString

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

    -- | GFF column 9. Feature attributes (see 'Attribute')
    , gff_attr :: ![Attribute]
  }
  deriving(Eq,Ord,Show)

-- | Construct a GffEntry from the full data of a parsed GFF line
gffEntry
  :: ByteString   -- ^ seqid
  -> ByteString   -- ^ source
  -> IntervalType -- ^ type
  -> Integer      -- ^ start
  -> Integer      -- ^ stop
  -> ByteString   -- ^ score
  -> Maybe Strand -- ^ strand
  -> ByteString   -- ^ phase
  -> [Attribute]  -- ^ attributes
  -> GffEntry
gffEntry seqid _ ftype start stop _ strand _ attr = 
  GffEntry {
      gff_seqid    = toShort seqid
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
  } = unsplit '\t'
    [
        fromShort seqid
      , "."
      , bshow ftype
      , bshow start
      , bshow stop
      , "."
      , maybe "." bshow strand
      , "."
      , (unsplit ';' . map bshow $ attr)
    ]


type GParser a = ByteString -> ReportS a

readInt :: GParser Integer
readInt s = case {-# SCC "gffEntry_integer" #-} readInteger s of
  Just (x,"") -> pass' x
  _ -> fail' $ "GffParse: expected integer, found '" ++ s ++ "'"

readType :: GParser IntervalType
readType s = {-# SCC "gffEntry_type" #-} case s of
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
  x                 -> pass' $ Other $ toShort x

readStrand :: GParser (Maybe Strand)
readStrand s = {-# SCC "gffEntry_strand" #-} case s of
  "+" -> pass' $ Just Plus
  "-" -> pass' $ Just Minus
  "." -> pass' Nothing -- The distinction between '.' and '?' is a nuance I
  "?" -> pass' Nothing -- will ignore for now
  v   -> fail' $ "GffParse: expected strand from set [+-.?], found '" ++ v

readAttributes :: GParser [Attribute]
readAttributes s =
  ( {-# SCC "readAttributes_top" #-}  sequence
    $ map ({-# SCC "readAttributes_toPair" #-} toPair)
    $ map ({-# SCC "readAttributes_mapsplitEq" #-} split '=')
    $ {-# SCC "readAttributes_mapsplitSC" #-} split ';'
    $ s
  ) >>= warnIfTagsRepeat >>= {-# SCC "readAttributes_toAttr" #-}(sequence . map toAttr) where

  toPair :: [ByteString] -> ReportS (ByteString, ByteString)
  toPair []             = fail' $ "GffAttribute: expected attribute (<tag>=<val>), found ''"
  toPair [tag,val]      = pass' (tag , val)
  toPair [val]          = pass' (""  , val)
  toPair fs             = fail' $ concat 
    ["GffAttribute: expected attribute (<tag>=<val>), found '", unsplit '=' fs, "'"]

  toAttr :: (ByteString, ByteString) -> ReportS Attribute
  -- single value entries
  toAttr ("ID"           , v) = pass' $ AttrID          $ toShort v
  toAttr ("Name"         , v) = pass' $ AttrName        $ toShort v
  toAttr ("Target"       , v) = pass' $ AttrTarget      $ toShort v
  toAttr ("Gap"          , v) = pass' $ AttrGap         $ toShort v
  toAttr ("Derives_from" , v) = pass' $ AttrDerivesFrom $ toShort v
  -- multiple value entries
  toAttr ("Alias"         , vs) = pass' $ AttrAlias        $ map toShort $ split ',' vs
  toAttr ("Parent"        , vs) = pass' $ AttrParent       $ map toShort $ split ',' vs
  toAttr ("Note"          , vs) = pass' $ AttrNote         $ map toShort $ split ',' vs
  toAttr ("Dbxref"        , vs) = pass' $ AttrDbxref       $ map toShort $ split ',' vs
  toAttr ("Ontology_term" , vs) = pass' $ AttrOntologyTerm $ map toShort $ split ',' vs
  -- boolean entries
  toAttr ("Is_circular", "true"  ) = pass' $ AttrIsCircular True
  toAttr ("Is_circular", "false" ) = pass' $ AttrIsCircular False
  toAttr ("Is_circular", _       ) = fail' "GFFAttribute: Is_circular tag must be either 'true' or 'false'"
  -- untagged entry
  -- TODO interpret untagged value as an ID if no ID is provided The spec does
  -- not require this, but I should add it in to handle the shit _certain_
  -- programs throw at us.
  toAttr ("", v) = pass' $ AttrUntagged (toShort v)
  -- other entries
  -- TODO: Note if the tag is upper case, since these should be
  --       reserved for future use
  toAttr (v, t) = pass' $ AttrUserDefined (toShort v) (toShort t)

  warnIfTagsRepeat :: [(ByteString, ByteString)] -> ReportS [(ByteString, ByteString)]
  warnIfTagsRepeat ps = case mapMaybe headMay . map (drop 1) . group . sort . map fst $ ps of
    [] -> pass' ps
    es -> pass' ps >>= warn' ("GFFAttribute: each tag may appear only once, offending tag(s): [" ++ tags ++ "]") where
      tags = unsplit ',' es


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
                          
  map ( {-# SCC "readGFF_splitTAB" #-} split '\t') . -- Break tests by line and TAB. NOTE:

  
  {-# SCC "readGFF_lines" #-} split '\n'                     -- this allows space in fields

  where
    toGff (_, [c1,c2,c3,c4,c5,c6,c7,c8,c9])
      = {-# SCC "gffEntry" #-} gffEntry
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
