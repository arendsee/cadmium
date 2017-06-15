{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FunctionalDependencies #-}

module Fagin.Prelude (
    module CorePrelude
  , module Safe
  -- monoid operator
  , (++)
  -- semigroup (with <> operator)
  , Semigroup(..) 
  , BShow(..)
  , Sequence(..)
  , P.show
  , DL.filter
  , DL.zip
  , DL.sort
  , DL.group
  , DL.lookup
  , map
  , DL.concatMap
  -- Control
  , CM.sequence
  , CM.mapM
  , DF.fold
  , DF.foldr
  , DF.foldr'
  , DF.foldl
  , DF.foldl'
  -- ByteString IO
  , print
  , DBC.readFile
  , DBC.writeFile
  , DBC.appendFile
  , DBC.getContents
  , DBC.putStr
  , DBC.putStrLn
  , DBC.interact
  , DBC.hGet
  , DBC.hGetLine
  , DBC.hGetContents
  , DBC.hPut
  , DBC.hPutStr
  , DBC.hPutStrLn
  -- ByteString textual functions
  , DBC.split
  ,     unsplit
  , DBC.words
  , DBC.unwords
  , DBC.lines
  , DBC.unlines
) where

{-|

For Prelude, I will use the CorePrelude, which is a fairly minimal Prelude that
has most of what I want. The main change I add is to hide the Text functions
and import Data.BytesString.Char8. Effectively this allows ByteString to act as
the primary String class. 

'ByteString' uses half the memory of 'Text' at the cost of handling only ASCII
characters. This is a reasonable tradeoff here since memory is likely to be
limiting and the data are nearly guaranteed to be only ASCII.

I also hide 'fail' from standard Prelude. I do this because I redefine it (like
a fool) in Fagin.Report.

-}

import CorePrelude hiding (
    -- hide Text functions
      putStr   
    , putStrLn
    , getArgs
    , terror
    -- hide String functions
    , print
    -- hide unsafe functions
    , fail
    , (<>)
    )

import Safe (headMay)

import qualified Prelude as P
import qualified Data.ByteString.Char8 as DBC
import qualified Data.List as DL
import qualified Control.Monad as CM
import qualified Data.Foldable as DF

class Semigroup a where
  infixr 6 <>
  (<>) :: a -> a -> a

infixr 6 ++
(++) :: Monoid a => a -> a -> a
(++) = mappend

class (Monoid f) => Sequence f e | f -> e where
  length      :: f -> Int
  reverse     :: f -> f
  intersperse :: e -> f -> f
  transpose   :: [f] -> [f]
  replicate   :: Int -> e -> f
  takeWhile   :: (e -> Bool) -> f -> f
  dropWhile   :: (e -> Bool) -> f -> f
  concat      :: [f] -> f
  isPrefixOf  :: f -> f -> Bool
  isSuffixOf  :: f -> f -> Bool
  take        :: Int -> f -> f
  drop        :: Int -> f -> f
  uncons      :: f -> Maybe (e,f)

instance Sequence ByteString Char where
  length      = DBC.length
  reverse     = DBC.reverse
  intersperse = DBC.intersperse
  transpose   = DBC.transpose
  replicate   = DBC.replicate
  takeWhile   = DBC.takeWhile
  dropWhile   = DBC.dropWhile
  concat      = DBC.concat     
  isPrefixOf  = DBC.isPrefixOf
  isSuffixOf  = DBC.isSuffixOf
  take        = DBC.take
  drop        = DBC.drop
  uncons      = DBC.uncons

instance (Eq a) => Sequence [a] a where
  length      = DL.length
  reverse     = DL.reverse
  intersperse = DL.intersperse
  transpose   = DL.transpose
  replicate   = DL.replicate
  takeWhile   = DL.takeWhile
  dropWhile   = DL.dropWhile
  concat      = DL.concat     
  isPrefixOf  = DL.isPrefixOf
  isSuffixOf  = DL.isSuffixOf
  take        = DL.take
  drop        = DL.drop
  uncons      = DL.uncons

map :: (Functor f) => (a -> b) -> f a -> f b
map = P.fmap

unsplit :: Char -> [ByteString] -> ByteString
unsplit c = DBC.intercalate (DBC.singleton c)

-- Set default IO to ByteString

print :: (BShow a) => a -> IO ()
print = DBC.putStrLn . bshow 

class (Show a) => BShow a where
  bshow :: a -> ByteString
  bshow = DBC.pack . P.show

instance BShow Int

instance BShow Integer

instance BShow String where
  bshow = DBC.pack

instance BShow ByteString where
  bshow = id

instance BShow a => BShow [a] where
  bshow = DBC.unlines . map bshow
