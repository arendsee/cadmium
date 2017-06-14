{-# LANGUAGE TypeFamilies #-}

module Fagin.Prelude (
    module ClassyPrelude
  , Textual(..)
  , BShow(..)
) where

{-|

The Prelude defined here is based almost completely off 'classy-prelude'. The
one change is the switch from 'Text' to 'ByteString' as the main string type.
The 'Textual' typeclass in Classy is defined as sequences of 'Char', I relax
this restriction.

'ByteString' uses half the memory of 'Text' at the cost of handling only ASCII
characters. This is a reasonable tradeoff here since memory is likely to be
limiting and the data are nearly guaranteed to be only ASCII.

I also hide 'fail' from standard Prelude. I do this because I redefine it (like
a fool) in Fagin.Report.

-}

import qualified Data.ByteString.Char8 as BC
import qualified Data.Char as DC
import ClassyPrelude hiding(
                               Textual(..)
                             , Utf8
                             , fail -- redefined in Report
                           )

class (IsSequence t, IsString t) => Textual t where
  words      :: t -> [t]
  unwords    :: (Element seq ~ t, MonoFoldable seq) => seq -> t
  lines      :: t -> [t]
  unlines    :: (Element seq ~ t, MonoFoldable seq) => seq -> t
  toLower    :: t -> t
  toUpper    :: t -> t
  toCaseFold :: t -> t
  breakWord  :: t -> (t, t)
  breakLine  :: t -> (t, t)

instance Textual ByteString where
  words      = BC.words
  unwords    = BC.unwords . otoList
  lines      = BC.lines
  unlines    = BC.unlines . otoList
  toLower    = BC.map DC.toLower
  toUpper    = BC.map DC.toUpper
  toCaseFold = toLower
  breakWord  = fmap (BC.dropWhile DC.isSpace) . BC.break DC.isSpace
  breakLine  = fmap (BC.dropWhile isNewline) . BC.break isNewline
    where
    isNewline = \c -> c == '\n' || c == '\r'

class (Show a) => BShow a where
  bshow :: a -> ByteString
  bshow = BC.pack . show

instance BShow Int

instance BShow Integer

instance BShow String where
  bshow = BC.pack

instance BShow ByteString where
  bshow = id

instance BShow a => BShow [a] where
  bshow = unlines . map bshow
