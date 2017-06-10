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
) where

import qualified Data.Text as T
import           Data.Monoid ((<>))

import Fagin.Interval
import Fagin.Error

data IntervalType
  = MRna
  | CDS
  | Exon
  | Gene
  | Other T.Text
  deriving(Show,Eq,Ord)

data GffEntry = GffEntry {
      gff_seqid    :: T.Text
    , gff_type     :: IntervalType
    , gff_strand   :: Maybe Strand
    , gff_interval :: Interval
    , gff_attr     :: [(T.Text, T.Text)]
  }
  deriving(Show,Eq,Ord)

type Attribute = (T.Text, T.Text)


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

readAttribute :: GParser [Attribute]
readAttribute = sequence . map asAttr . map (T.splitOn "=") . T.splitOn ";" where
  asAttr :: [T.Text] -> ThrowsError Attribute
  asAttr []             = Left $ GffExpectAttribute ""
  asAttr [tag,val]      = Right (tag , val)
  asAttr [val]          = Right (""  , val)
  asAttr fs             = Left $ GffExpectAttribute $ T.intercalate "=" fs


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
      = case (readType typ, readInt a, readInt b, readStrand str, readAttribute attr) of
        -- If everything is good, make an entry
        (Right typ', Right a', Right b', Right str', Right attr') ->
          [ Right $
              GffEntry {
                  gff_seqid    = chr
                , gff_type     = typ'
                , gff_strand   = str'
                , gff_interval = Interval a' b'
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
