{-# LANGUAGE OverloadedStrings #-}

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
  ""     -> Left  GffNoType
  "mRNA" -> Right MRna
  "CDS"  -> Right CDS
  "exon" -> Right Exon
  "Gene" -> Right Gene
  x      -> Right $ Other x

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
