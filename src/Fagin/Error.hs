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
  show ModelInvalidParent = "Parent tag matches no entry Id" -- TODO parameterize
  show ModelExpectParent = "Entry of this type ought to have a parent" --TODO parameterize
  show ModelStrandMissing = "Gene model CDS and exon entries do not specify strand"
  show ModelStrandMismatch = "Gene model has contradictory strand specifications"
  -- GFF
  show (GffInvalidRowNumber xs)
    | (length xs) < 9 = "Too few columns"
    | (length xs) > 9 = "Too many columns"
    | otherwise       = "Well shucks, that shouldn't have happened"
  show (GffExpectInteger v)     = "Expected integer, found '" ++ T.unpack v ++ "'"
  show (GffExpectAttribute msg) = "Expected attribute (<tag>=<val>), found '" ++ T.unpack msg ++ "'"
  show (GffExpectStrand v)      = "Expected strand ('+', '-' or '.'), found '" ++ T.unpack v ++ "."
  show GffNoFeatures            = "No features found (empty file)"
  show GffNoType                = "Exected <type> in column 3, found nothing"
  show (MultiError [])       = ""
  show (MultiError es)       = " - " ++ (L.intercalate "\n - " . map show $ es)
  show NoError               = ""
  -- contectualize a GFF error with line information
  show (GffLineError i (MultiError es)) =
    "line " ++ show i ++ ":\n - " ++
    concatMap (\s -> " - " ++ show s ++ "\n") es
  show (GffLineError _ NoError) = ""
  show (GffLineError i e) = "line " ++ show i ++ ": " ++ show e
