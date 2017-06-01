module Fagin.Interval (
    Strand(..)
  , Interval(..)
) where

data Strand
  = Plus
  | Minus
  deriving(Show)

data Interval = Interval Integer Integer (Maybe Strand) deriving(Show)
