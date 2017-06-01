{-# LANGUAGE OverloadedStrings #-}

module Fagin.Gff (
    readGff
  , IntervalType(..)
  , Attribute(..)
  , GffEntry(..)
) where

import qualified Data.Attoparsec.Text as AT
import qualified Data.Text as T

import Fagin.Interval

data IntervalType
  = MRna
  | CDS
  | Exon
  | Gene
  deriving(Show)

data Attribute
  = Parent   T.Text
  | Id       T.Text
  | Name     T.Text
  | Untagged T.Text
  | Other    T.Text T.Text
  | Nil
  deriving (Show)

data GffEntry = GffEntry {
      gff_seqid    :: T.Text
    , gff_type     :: Maybe IntervalType
    , gff_interval :: Interval
    , gff_parent   :: Maybe T.Text
    , gff_id       :: Maybe T.Text
    , gff_name     :: Maybe T.Text
    , gff_meta     :: [(T.Text, T.Text)]
  }
  deriving (Show)

parseField :: AT.Parser T.Text
parseField = AT.takeWhile (\c -> c /= '\t')

parseComment :: AT.Parser T.Text
parseComment = do
  _ <- AT.char '#'
  s <- AT.takeWhile (\c -> c /= '\n')
  _ <- AT.char '\n'
  return s

parseSeqid :: AT.Parser T.Text
parseSeqid = do
  s <- parseField
  return s

parseType :: AT.Parser (Maybe IntervalType)
parseType = do
  s <- AT.choice
    [
      AT.string "mRNA"
    , AT.string "CDS"
    , AT.string "exon"
    , AT.string "gene"
    ]
  return $ case s of
    "mRNA" -> Just MRna
    "CDS"  -> Just CDS
    "exon" -> Just Exon
    "gene" -> Just Gene
    _      -> Nothing

parseStrand :: AT.Parser (Maybe Strand)
parseStrand = do
  s <- AT.anyChar
  return $ case s of
    '-' -> Just Minus
    '+' -> Just Plus
    _   -> Nothing

parseAttributePair :: AT.Parser Attribute
parseAttributePair = do
  a <- AT.takeWhile (\c -> c /= '=')
  _ <- AT.char '='
  b <- AT.takeWhile (\c -> c /= ';')
  return $ case (a,b) of
    ("Parent" , val) -> Parent val
    ("Id"     , val) -> Id     val
    ("Name"   , val) -> Name   val
    (tag      , val) -> Other  tag val

parseAttributeUntagged :: AT.Parser Attribute
parseAttributeUntagged = do
  x <- AT.takeWhile (\c -> c /= ';')
  return $ case x of
    "" -> Nil
    s  -> Untagged s

parseAttribute :: AT.Parser Attribute
parseAttribute = AT.choice [parseAttributePair, parseAttributeUntagged]

parseMeta :: AT.Parser [Attribute]
parseMeta = AT.sepBy parseAttribute (AT.char ';')

parseGffEntry :: AT.Parser GffEntry
parseGffEntry = do
  chr <- parseSeqid
  _   <- AT.char '\t' 

  _   <- parseField
  _   <- AT.char '\t' 

  typ <- parseType
  _   <- AT.char '\t' 

  a   <- AT.decimal
  _   <- AT.char '\t' 

  b   <- AT.decimal
  _   <- AT.char '\t' 

  _   <- parseField
  _   <- AT.char '\t'

  str <- parseStrand
  _   <- AT.char '\t'

  _   <- AT.anyChar
  _   <- AT.char '\t'

  att <- parseMeta

  _   <- AT.many' AT.space

  return GffEntry {
      gff_seqid    = chr
    , gff_type     = typ
    , gff_interval = Interval a b str
    , gff_parent   = findParent att
    , gff_id       = findId att
    , gff_name     = findName att
    , gff_meta     = findMeta att
  } where

    findParent :: [Attribute] -> Maybe T.Text
    findParent []            = Nothing
    findParent (Parent s:_)  = Just s
    findParent (_:xs)        = findParent xs

    findId :: [Attribute] -> Maybe T.Text
    findId []       = Nothing
    findId (Id s:_) = Just s
    findId (_:xs)   = findId xs

    findName :: [Attribute] -> Maybe T.Text
    findName []         = Nothing
    findName (Name s:_) = Just s
    findName (_:xs)     = findName xs

    findMeta :: [Attribute] -> [(T.Text,T.Text)]
    findMeta [] = []
    findMeta (Other t v:xs) = [(t,v)] ++ findMeta xs
    findMeta _ = []

parseGff :: AT.Parser [GffEntry]
parseGff = do
  _   <- AT.many' parseComment
  gff <- AT.many1 parseGffEntry
  return gff

readGff :: T.Text -> Either String [GffEntry]
readGff = AT.parseOnly parseGff
