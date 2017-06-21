import Fagin.Prelude
import Fagin.Report
import Fagin.Interval

data Interval = Interval Start Stop

newtype ScaffoldID = ScaffoldID ByteString
newtype GeneID     = GeneID     ByteString
newtype ModelID    = ModelID    ByteString
newtype ParentID   = ParentID   ByteString

newtype DerivedFrom = DerivedFrom ByteString

data Genome    = Genome    [Scaffold]
data Scaffold  = Scaffold  ScaffoldID [Gene]
data Gene      = Gene      GeneID [GeneModel]
data GeneModel = GeneModel ModelID Orf Mrna Strand
data Orf       = Orf       [Interval]
data Trans     = Mrna      MRnaID [Interval]

data GeneEntry = GeneEntry GeneID
data MrnaEntry = MrnaEntry ModelID  Interval Strand [ParentID]
data ExonEntry = ExonEntry          Interval Strand [ParentID]
data CDSEntry  = CDSEntry           Interval Strand [ParentID] [DerivedFrom]


-- G0 - Initial graph
--
--    e_ ---n,mt--> m_ --m,p'--> g_
--                  ^            ^
--                 /            /
--    c_ ---cp,m--'----1,1-----'
--
--    where
--    e_ := GffEntry     n  := number of exons in a transcript
--    m_ := GffMRna      m  := number of transcripts in a model
--    c_ := GffCDS       t  := number of mRNA spliced into a transcript (transplicing)
--    g_ := GffGene      c  := number of CDS in a ORF
--                       p  := number of ORFs in a transcript
--                       p' := number of ORFs in a model


-- G1 - collapse exons into transcripts (assume no trans-splicing)
--
--            T ---m,p'--> g_
--            ^            ^
--         cp,m           /
--          /            /
--    c_ --'----1,1-----'
collapseExons ::  [ExonEntry] -> MrnaEntry -> Transcript

-- G2 - collapse CDS into ORFs
--
--                  T ---m,p'--> g_
--                  ^            ^
--               1,p            /
--               /             /
--        C ----'-----1,1-----'
collapseCDS :: [CDSEntry] -> Transcript -> ORF


-- G3 - infer protein names
--
--        C --n,1--> T ---m,p'--> g_
inferNames :: ORF -> Transcript ->


-- toMRna :: [MRnaEntry] -> [ExonEntry] -> ReportS [MRna]
-- toMRna = undefined
--
-- toGffRow :: [ByteString] -> ReportS GffRow
-- toGffRow = undefined


-- partitionTypes :: [GffRow] -> ([GeneEntry], [MRnaEntry], [ExonEntry], [CDSEntry])
-- partitionTypes = undefined


-- toMatrix :: ByteString -> ReportS [GffRow]
-- toMatrix =
--   sequence
--   map toGffRow                 . -- to record, catch row number errors
--   map (split '\t')             . -- Break tests by line and TAB.
--   filter (\s -> length s == 0) . -- Remove empty lines
--   filter (isPrefixOf "#")      . -- Remove comments
--   split '\n'                     -- this allows space in fields


-- readGff :: ByteString -> ReportS [Gene]
-- readGff = map toFeature . toMatrix




-- data GffRow = GffRow {
--       gff_seqid  :: ByteString
--     , gff_source :: ByteString
--     , gff_type   :: ByteString
--     , gff_start  :: ByteString
--     , gff_stop   :: ByteString
--     , gff_score  :: ByteString
--     , gff_strand :: ByteString
--     , gff_phase  :: ByteString
--     , gff_attr   :: ByteString
--   }



---- Ultimate Goal
-- Genome = Genome [(Seqid, [Gene])]
-- Gene = Gene [(ModelId, [GeneModel])]
-- GeneModel = GeneModel {
--     getTrs :: [Exon]
--   , getTrl :: [CDS]
--   , getStr :: Strand
-- }


--   text -> Genome [(Seqid, [Gene])]
--  ==================================================================
--   text
-- ------------------------------------------------------------------
--   [(seqid,...)]
-- ------------------------------------------------------------------
--   [(seqid, [(type,interval,strand,attr)])]
--   ==> Geneome [(Seqid, f x)]
-- ------------------------------------------------------------------
--   f :: [(type,interval,strand,attr)] -> [Gene]
--   
-- ------------------------------------------------------------------
--
--     [Mrna] [CDS] [Exon] [Gene]
--       where
--         Mrna :: ID [Parent]




-- type ID = ConstantString
--
-- data Edge = PartOf Node | DerivedFrom Node | PartOf' ID | DerivedFrom' ID
--
-- data Node = Node Feature [Edge]
--
-- data Feature
--   = Gene GeneId [Attribute]
--   | Mrna Start Stop Strand SeqId [Attribute]
--   | Exon Start Stop Strand
--   | CDS  Start Stop Strand

-- data GeneMode = GeneModel {
--       model_mrna :: MrnaId
--     , model_cds :: [Interval]
--   }
--
-- data GeneModel = GeneModel {
--       model_seqid  :: !SeqId
--     , model_geneid :: !(Maybe GeneId) -- based either on CDS:grandparent or CDS:Derived_from
--     , model_mrnaid :: !(Maybe MrnaId) -- based on CDS:Parent and Exon:Parent
--     , model_cds    :: ![Interval]
--     , model_exon   :: ![Interval]
--     , model_strand :: !Strand
--   }
--
-- type ChildId  = ConstantString
-- type ParentId = ConstantString
-- type GeneId   = ConstantString
-- type MrnaId   = ConstantString
-- type SeqId    = ConstantString
--
-- type Start = Int
-- type Stop  = Int
--
-- data Strand = Plus | Minus | Unset
--
-- newtype Gene = Gene GeneId [Attribute]
-- newtype Mrna = Mrna Start Stop Strand (Maybe GeneId) SeqId [Attribute]
-- newtype Exon = Exon Start Stop Strand MrnaId
-- newtype CDS  = CDS  Start Stop Strand MrnaId (Maybe DerivedFrom)
--
-- data Transcript = Transcript [Interval] Strand SeqId [Parent] [Attribute]
--
-- data Translant  = Translant [Interval] MrnaID (Maybe DerivedID)
--
-- partitionTypes :: [[ByteString]] -> ([Gene],[MRna],[Exon],[CDS])
-- partitionTypes = undefined
--
-- -- | Some exons will have multiple parents, I unravel this particular demon
-- -- here. This step permanently removes Exon types from consideration, wrapping
-- -- them into transcipts.
-- toTranscript :: [MRna] -> [Exon] -> HashMap MrnaId Transcript
--
-- toTranslant :: [CDS] -> [Translant]
--
-- toHash :: [[ByteString]] -> HashMap ParentId

