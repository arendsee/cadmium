# File -> AAStringSet
load_fasta_aa <- function(filename){
  require(Biostrings)
  readAAStringSet(filename)
}

# File -> DNAStringSet
load_fasta_fna <- function(filename){
  require(Biostrings)
  readDNAStringSet(filename)
}

# File -> GFF_raw
load_gff_table <- function(filename){ }

# File -> [Char]
load_character_vector <- function(filename){ }

# File -> GI
load_nstring <- function(filename){ }

# File -> Seqinfo -> Seqinfo -> (SI, unassembled)
load_search_intervals <- function(filename){
  synder::load_si(filename)
}

# File -> Seqinfo -> Seqinfo -> (GI, GI)
load_synteny_map <- function(filename){
  synmap <- synder::load_synmap(filename)
}

# File -> (Species, Seqid, Length)
load_scaffold_lengths <- function(filename){ }

# File -> Tree
load_tree <- function(treefile){
  ape::read.tree(treefile)
}

# Named -> Named
reduce_fasta_headers <- function(fasta){
  # Remove any comments that follow the sequence name in the header
  names(fasta) <- gsub(' .*', '', names(fasta))
  fasta
}

# (Species, Seqid, Length) -> Seqinfo
get_seqinfo <- function(scaflen){ }

# GFF_raw -> Seqinfo -> GI
get_annotated_gi <- function(gff, seqinfo){ }

# Position -> Position -> Seqid -> GI
make_GI <- function(starts, stops, seqids){ }

# GI -> String -> GI
gi_add_strand <- function(gi, strand){ }

# GI -> Seqinfo -> GI
gi_add_seqinfo <- function(gi, seqinfo){ }

# GI -> Table -> GI
gi_add_meta <- function(gi, meta){
  mcols(gi) <- meta
  gi
}
