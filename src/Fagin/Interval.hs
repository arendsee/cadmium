module Fagin.Interval (
    Strand(..)
  , Interval(..)
) where

data Strand
  = Plus
  | Minus
  deriving(Show,Ord,Eq)

data Interval = Interval Integer Integer deriving(Show,Ord,Eq)
