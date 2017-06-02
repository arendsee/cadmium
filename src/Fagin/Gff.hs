{-# LANGUAGE OverloadedStrings #-}

module Fagin.Gff (
    readGff
  , IntervalType(..)
  , Attribute(..)
  , GffEntry(..)
  , CGffError
) where

import qualified Data.Text as T
import qualified Data.List as L
import qualified Data.Either as E
import           Data.Monoid ((<>))

import Fagin.Interval

data IntervalType
  = MRna
  | CDS
  | Exon
  | Gene
  | Other T.Text
  deriving(Show)

data Attribute
  = Parent   T.Text
  | Id       T.Text
  | Name     T.Text
  | Untagged T.Text
  | Meta     T.Text T.Text
  | Nil
  deriving (Show)

data GffEntry = GffEntry {
      gff_seqid    :: T.Text
    , gff_type     :: IntervalType
    , gff_interval :: Interval
    , gff_parent   :: Maybe T.Text
    , gff_id       :: Maybe T.Text
    , gff_name     :: Maybe T.Text
    , gff_meta     :: [Attribute]
  }
  deriving (Show)

type Row = [T.Text]
type Message = T.Text

data GffError
  = NoError
  | NoType
  | NoFeatures
  | InvalidRowNumber Row
  | ExpectInteger T.Text
  | ExpectAttribute Message
  | ExpectStrand T.Text
  | MultiError [GffError]

instance Monoid GffError where
  mempty = NoError
  mappend NoError NoError = NoError
  mappend NoError y       = y
  mappend x       NoError = x
  mappend (MultiError xs) (MultiError ys) = MultiError (xs  ++ ys)
  mappend (MultiError xs) x               = MultiError (xs  ++ [x])
  mappend x               (MultiError ys) = MultiError ([x] ++ ys)
  mappend x               y               = MultiError [x,y]

instance Show GffError where
  show (InvalidRowNumber xs)
    | (length xs) < 9 = "Too few columns"
    | (length xs) > 9 = "Too many columns"
    | otherwise       = "Well shucks, that shouldn't have happened"
  show (ExpectInteger v)     = "Expected integer, found '" ++ T.unpack v ++ "'"
  show (ExpectAttribute msg) = "Expected attribute (<tag>=<val>), found '" ++ T.unpack msg ++ "'"
  show (ExpectStrand v)      = "Expected strand ('+', '-' or '.'), found '" ++ T.unpack v ++ "."
  show NoFeatures            = "No features found (empty file)"
  show NoType                = "Exected <type> in column 3, found nothing"
  show (MultiError [])       = ""
  show (MultiError es)       = " - " ++ (L.intercalate "\n - " . map show $ es)
  show NoError               = ""


data CGffError = CGffError Integer GffError

instance Show CGffError where
  show (CGffError i (MultiError es)) =
    "line " ++ show i ++ ":\n - " ++
    concatMap (\s -> " - " ++ show s ++ "\n") es
  show (CGffError _ NoError) = ""
  show (CGffError i e) = "line " ++ show i ++ ": " ++ show e

reerror :: [Either l r] -> Either l [r]
reerror es = case E.partitionEithers es of
  ([],rs)   -> Right rs
  ((l:_),_) -> Left l


type GParser a = T.Text -> Either GffError a

readInt :: GParser Integer
readInt s = case reads (T.unpack s) :: [(Integer,String)] of
  [(x,"")] -> Right x
  _        -> Left $ ExpectInteger s

readType :: GParser IntervalType
readType s = case s of
  ""     -> Left  NoType
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
  v   -> Left $ ExpectStrand v

readAttribute :: GParser [Attribute]
readAttribute = reerror . map asAttr . map (T.splitOn "=") . T.splitOn ";" where
  asAttr :: [T.Text] -> Either GffError Attribute
  asAttr []             = Right Nil
  asAttr ["Parent",val] = Right $ Parent val
  asAttr ["ID",val]     = Right $ Id val
  asAttr ["Name",val]   = Right $ Name val
  asAttr [tag,val]      = Right $ Meta tag val
  asAttr [val]          = Right $ Untagged val
  asAttr fs             = Left $ ExpectAttribute $ T.intercalate "=" fs


type GIFilter = (Integer,[T.Text]) -> Bool

comment :: GIFilter
comment (_,(x:_)) = T.isPrefixOf (T.singleton '#') x
comment _ = True

empty :: GIFilter
empty (_,[]) = True
empty _      = False

readGff :: T.Text -> Either CGffError [GffEntry]
readGff =
  reerror                . -- Merge errors, die on first failure

  toGff                  . -- Parse GFF entry and report errors

  filter (not . empty)   . -- Filter out lines that are either empty
  filter (not . comment) . -- or start with a comment (#) character

  zip [1..]              . -- Add line numbers. This must precede filters
                           -- so line numbering in error messages is correct

  map (T.splitOn "\t")   . -- Break tests by line and TAB. NOTE:
  T.lines                  -- this allows space in fields
  where
    toGff :: [(Integer,[T.Text])] -> [Either CGffError GffEntry]
    toGff ((i, [chr, _, typ, a, b, _, str, _, attr]):xs)
      = case (readType typ, readInt a, readInt b, readStrand str, readAttribute attr) of
        (Right typ', Right a', Right b', Right str', Right attr') ->
          [ Right $
              GffEntry {
                  gff_seqid    = chr
                , gff_type     = typ'
                , gff_interval = Interval a' b' str'
                , gff_parent   = Nothing
                , gff_id       = Nothing
                , gff_name     = Nothing
                , gff_meta     = attr'
              }
          ] ++ toGff xs
        (typ', a', b', str', attr') -> [ Left (CGffError i err) ] where
          l :: Either GffError a -> GffError
          l (Right _) = mempty
          l (Left  e) = e
          err = l typ' <> l a' <> l b' <> l str' <> l attr'
    toGff [] = []
    toGff ((i,fs):_) = [ Left $ CGffError i $ InvalidRowNumber fs ]
