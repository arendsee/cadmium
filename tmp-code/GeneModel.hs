
data GeneModel = Model {
  model_parent :: Maybe Id
  model_id     :: Id
  model_cds    :: [Pos]
  model_exon   :: [Pos]
  model_strand :: Maybe Strand
}

