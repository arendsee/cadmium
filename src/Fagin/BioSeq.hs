{-# LANGUAGE OverloadedStrings #-}

module Fagin.BioSeq (BioSeq(..)) where

import qualified Data.ByteString.Char8 as B

-- import Fagin.Report
--
-- data DNA
-- data AA
-- data Unknown

type RawSeq = B.ByteString
type SeqId = B.ByteString
type Desc = B.ByteString


data BioSeq a
  = BioSeq {
      seqid   :: SeqId
    , seqdesc :: Maybe Desc
    , seqstr  :: RawSeq
  }

-- instance Show (BioSeq a) where
--   show = show . toFasta 
-- 
-- seqslice :: BioSeq a -> Int -> Int -> ThrowsError (BioSeq a)
-- seqslice i j b
--   | len == j  = Right $ BioSeq newseqid Nothing newseqstr
--   | len == 0  = Left $ IndexError (msg ++ ", i > seqlen")
--   | otherwise = Left $ IndexError (msg ++ ", j > seqlen")
--   where
--     newseqstr = B.take j . snd . B.splitAt i . seqstr $ b
--     newseqid = seqid b -- ++ (B.pack "_[") ++ (B.pack . show) i ++ (B.singleton ',') ++ (B.pack . show) j ++ (B.singleton ']')
--     len = length s
--     msg = "cannot extract slice (i=" ++
--           show i ++ ", j=" ++ show j ++ ")" ++
--           " from '" ++ seqid b ++ "'"
-- 
-- fromFasta :: B.ByteString -> [BioSeq Unknown]
-- fromFasta = map B.split '\n' . B.split '>'
-- 
-- toFasta :: BioSeq a -> B.ByteString
-- toFasta s = concat [">", id', " ", desc', "\n", seq'] where
--   id'   = seqid s
--   desc' = fromMaybe " " (seqdesc s)
--   seq'  = seqstr s
