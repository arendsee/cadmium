module Fagin (fullData) where

import Fagin.Gff (openGff)

fullData :: String -> String
fullData _ = openGff "Might I have some more?"
