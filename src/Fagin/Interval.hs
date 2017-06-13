module Fagin.Interval (
    Strand(..)
  , Interval(..)
  , intervalSpan
) where

data Strand
  = Plus
  | Minus
  deriving(Ord,Eq)

instance Show Strand where
  show Plus  = "+"
  show Minus = "-"

data Interval = Interval Integer Integer deriving(Show,Ord,Eq)

intervalSpan :: Interval -> Interval -> Interval
intervalSpan (Interval a1 b1) (Interval a2 b2) = Interval (min a1 a2) (max b1 b2) 
