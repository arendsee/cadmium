module Fagin (setup) where

data Id
  = Scaf  String
  | Model String
  | Gene  String
  | MRna  String
  | Taxon Integer

data Seq
  = Dna String
  | Pro String

data IntervalType = MRna | CDS | Exon | Gene

data Strand = Plus | Minus

data Attribute
  = Parent   String
  | Id       String
  | Name     String
  | Untagged String
  | Other    String

data GffEntry = GffEntry {
  gff_type     :: IntervalType
  gff_interval :: Interval
  gff_parent   :: Maybe String
  gff_id       :: Maybe String
  gff_name     :: Maybe String
  gff_meta     :: [(String, String)]
}

data Gff = Gff [GffEntry]

type Pos = (Integer, Integer)

type Interval = (Id, Pos, Maybe Strand)

type Genome = [(Scaffold, DnaSeq)]

type Scaffold = Scaffold ScafId [GeneModel]

data GeneModel = Model {
  model_parent :: Maybe Id
  model_id     :: Id
  model_cds    :: [Pos]
  model_exon   :: [Pos]
  model_strand :: Maybe Strand
}

data FaginInput = FaginInput Focal [Target]
