module Interval (
    Strand(..)
  , Pos
  , Interval
) where

data Strand = Plus | Minus

type Pos = (Integer, Integer)

type Interval = (Id, Pos, Maybe Strand)
