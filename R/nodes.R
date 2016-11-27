h_load_genome   <- hwell('Genome')
h_load_synmap   <- hwell('Synmap')
h_species_tree  <- hwell('Tree')
h_load_gff      <- hwell('GFF')

h_get_seq_lengths <- hpipe('Genome -> SeqLengths')
h_annotate_gff    <- hpipe('SeqLengths -> GFF')

h_get_CDS         <- hpipe('Genome -> GFF -> CDS')
h_get_proteome    <- hpipe('CDS -> Pro_Seq')
h_get_transcripts <- hpipe('Genome -> GFF -> DNA_Seq')

h_get_genome_orfs_gff     <- hpipe('Genome -> GFF')
h_get_genome_orfs_dna_seq <- hpipe('GFF -> Genome -> CDS')
h_get_genome_orfs_pro_seq <- hpipe('CDS -> Pro_Seq')

h_get_seq <- hpipe('Genome -> GFF -> DNA_Seq')

h_subset_gff <- hpipe('GFF -> GFF')
h_subset_synmap <- hpipe('Synmap -> Synmap')


# --- Features

h_find_gen <- hnode('Pro_Seq -> Pro_Seq -> Bool')
h_find_cds <- hnode('CDS -> Pro_Seq -> Bool')
h_find_ind <- hnode('GFF -> GFF -> Bool')
h_find_nst <- hnode('GFF -> GFF -> Bool')
h_find_nuc <- hnode('DNA_Seq -> DNA_Seq -> Bool')
h_find_orf <- hnode('Pro_Seq -> Pro_Seq -> Bool')
h_find_res <- hnode('Synmap -> GFF -> Bool')
h_find_rna <- hnode('GFF -> GFF -> Bool')
h_find_scr <- hnode('Synmap -> Bool')
h_find_tec <- hnode('Synmap -> Bool')
h_find_trn <- hnode('Synmap -> Bool')
h_find_una <- hnode('Synmap -> Bool')
