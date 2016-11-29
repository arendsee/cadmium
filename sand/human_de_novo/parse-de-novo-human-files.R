# require(devtools)
# install_github('arendsee/sandr')

library(sandr)
library(tidyr)
library(dplyr)
library(magrittr)
library(readr)

gff <- sandr::read_sand('.', data_has_header=FALSE) %>%
  { .$orf_id <- paste0('orf_', 1:nrow(.)); . } %>%
  dplyr::mutate(
    exon_lengths = strsplit(exon_lengths, ','),
    exon_starts = strsplit(exon_starts, ',')
  ) %>%
  tidyr::unnest() %>%
  dplyr::mutate(
    phase = '.',
    exon_stops = as.numeric(exon_starts) + as.numeric(exon_lengths) - 1,
    attr = paste0(
      'transcript_start=', transcript_start, ';',
      'transcript_end=',   transcript_end,   ';',
      'CDS_start=',        CDS_start,        ';',
      'CDS_end=',          CDS_end,          ';',
      'label=',            label,            ';',
      'nexons=',           nexons,           ';',
      'gene_name=',        gene_name,        ';',
      'biotype=',          biotype,          ';',
      'aa_seq=',           aa_seq
    )
  ) %>%
  dplyr::select(
    chr=chromosome,
    orf_id,
    type,
    start=exon_starts,
    stop=exon_stops,
    phase,
    score,
    attr
  )

faa <- dplyr::select(gff, orf_id, attr) %>%
  dplyr::mutate(seq = sub('.*aa_seq=([^;]+).*', '\\1', attr)) %>%
  with(paste0('>', orf_id, '\n', seq, collapse='\n'))

readr::write_tsv(gff, 'human_de_novo.gff3')
readr::write_file(faa, 'human_de_novo.faa')
