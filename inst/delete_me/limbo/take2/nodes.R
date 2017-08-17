# h_load_genome   <- hwell('Genome')
# h_load_proteome <- hwell('Proteome')
# h_load_gff      <- hwell('GFF')
# h_species_tree  <- hwell('Tree')
#
# h_load_synmap   <- hpipe('SeqLengths -> SeqLengths -> Synmap')
# h_get_seq_lengths <- hpipe('Genome -> SeqLengths')
# h_annotate_gff    <- hpipe('SeqLengths -> GFF')
#
# h_get_CDS         <- hpipe('Genome -> GFF -> CDS')
# h_get_proteome    <- hpipe('CDS -> Pro_Seq')
# h_get_transcripts <- hpipe('Genome -> GFF -> DNA_Seq')
#
# h_get_genome_orfs_gff     <- hpipe('Genome -> GFF')
# h_get_genome_orfs_dna_seq <- hpipe('GFF -> Genome -> CDS')
# h_get_genome_orfs_pro_seq <- hpipe('CDS -> Pro_Seq')
#
# h_get_seq <- hpipe('Genome -> GFF -> DNA_Seq')
#
# h_subset_gff <- hpipe('GFF -> GFF')
# h_subset_synmap <- hpipe('Synmap -> Synmap')
#
# h_find_gen <- hpipe('Pro_Seq -> Pro_Seq -> Bool')
# h_find_cds <- hpipe('CDS -> Pro_Seq -> Bool')
# h_find_ind <- hpipe('GFF -> GFF -> Bool')
# h_find_nst <- hpipe('GFF -> GFF -> Bool')
# h_find_nuc <- hpipe('DNA_Seq -> DNA_Seq -> Bool')
# h_find_orf <- hpipe('Pro_Seq -> Pro_Seq -> Bool')
# h_find_res <- hpipe('Synmap -> GFF -> Bool')
# h_find_rna <- hpipe('GFF -> GFF -> Bool')
# h_find_scr <- hpipe('Synmap -> Bool')
# h_find_tec <- hpipe('Synmap -> Bool')
# h_find_trn <- hpipe('Synmap -> Bool')
# h_find_una <- hpipe('Synmap -> Bool')
