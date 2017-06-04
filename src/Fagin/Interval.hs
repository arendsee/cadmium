module Fagin.Interval (
    Strand(..)
  , Interval(..)
  , Pos
) where

data Strand
  = Plus
  | Minus
  deriving(Show,Ord,Eq)

type Pos = (Integer, Integer)

data Interval
  = Interval Pos (Maybe Strand)
  deriving(Show,Ord,Eq)
