module Fagin.Gff (openGff) where

import qualified Data.ByteString as B
import qualified Data.Attoparsec.ByteString as A

import Fagin.Interval

data IntervalType = MRna | CDS | Exon | Gene

data Attribute
  = Parent   String
  | Id       String
  | Name     String
  | Untagged String
  | Other    String

data GffEntry = GffEntry {
    gff_type     :: IntervalType
  , gff_interval :: Interval
  , gff_parent   :: Maybe String
  , gff_id       :: Maybe String
  , gff_name     :: Maybe String
  , gff_meta     :: [(String, String)]
}

data Gff = Gff [GffEntry]

openGff :: String -> String
openGff s = s
