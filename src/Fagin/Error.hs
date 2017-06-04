module Fagin.Error
(
    FaginError(..)
  , ThrowsError
) where

import qualified Data.Text as T
import qualified Data.List as L

type ThrowsError = Either FaginError
type LineNumber = Integer

data FaginError
  = NoError
  | MultiError [FaginError]
  -- Errors in gene models
  | ModelInvalidParent
  | ModelExpectParent
  | ModelStrandMissing
  | ModelStrandMismatch
  -- Errors in GFF files
  | GffNoType
  | GffNoFeatures
  | GffInvalidRowNumber [T.Text]
  | GffExpectInteger    T.Text
  | GffExpectAttribute  T.Text
  | GffExpectStrand     T.Text
  | GffLineError LineNumber FaginError
  deriving(Eq,Ord)

instance Monoid FaginError where
  mempty = NoError
  mappend NoError NoError = NoError
  mappend NoError y       = y
  mappend x       NoError = x
  mappend (MultiError xs) (MultiError ys) = MultiError (xs  ++ ys)
  mappend (MultiError xs) x               = MultiError (xs  ++ [x])
  mappend x               (MultiError ys) = MultiError ([x] ++ ys)
  mappend x               y               = MultiError [x,y]

instance Show FaginError where
  -- GeneModel
  show ModelInvalidParent  = "ModelInvalidParent: 'Parent' tag matches no entry 'Id' tag"
  show ModelExpectParent   = "ModelExpectParent: entry of this type ought to have a parent"
  show ModelStrandMissing  = "ModelStrandMissing: CDS and exon entries specify no strand"
  show ModelStrandMismatch = "ModelStrandMismatch: all elements of a gene model must be on the same strand"
  -- GFF
  show (GffInvalidRowNumber xs)
    | (length xs) < 9 = "GffInvalidRowNumber: too few columns"
    | (length xs) > 9 = "GffInvalidRowNumber: too many columns"
    | otherwise       = "GffInvalidRowNumber: well shucks, that shouldn't have happened"
  show (GffExpectInteger v)     = "GffExpectInteger: expected integer, found '" ++ T.unpack v ++ "'"
  show (GffExpectAttribute msg) = "GffExpectAttribute: expected attribute (<tag>=<val>), found '" ++ T.unpack msg ++ "'"
  show (GffExpectStrand v)      = "GffExpectStrand: expected strand ('+', '-' or '.'), found '" ++ T.unpack v ++ "."
  show GffNoFeatures            = "GffNoFeatures: no features found (empty file)"
  show GffNoType                = "GffNoType: exected <type> in column 3, found nothing"
  show (MultiError [])       = ""
  show (MultiError es)       = "MultiError:\n - " ++ (L.intercalate "\n - " . map show $ es)
  show NoError               = ""
  -- contectualize a GFF error with line information
  show (GffLineError i e) = "GFF line " ++ show i ++ ": " ++ show e
