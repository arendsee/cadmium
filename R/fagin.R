#' @importFrom methods isClass new
#' @importFrom graphics plot
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
#' @importFrom rmonad "%>>%" "%v>%" "%*>%" "%__%" "%v__%" "%||%" "%|>%" "%>_%"
#' @importFrom utils head tail
utils::globalVariables(c("%>%", ".", "%>>%", "%v>%", "%*>%", "%__%", "%v__%", "%||%", "%|>%", "%>_%"))
NULL

#' fagin: Trace the origins of orphan genes
#'
#' This documentation is a work in progress ...
#'
#' @section GFF input:
#'
#' Absolute path to a directory containing a GFF file for each species used in
#' the pipeline. This GFF file must contain at minimum mRNA and coding sequence
#' (CDS) features. All start and stop positions must be relative to the
#' reference genomes in FNA_DIR (see argument -n).
#'
#' The following must be true of all GFF files:
#' \itemize{
#'   \item Contain CDS, exon and mRNA entries
#'   \item All CDS and exon fields contain Parent tags in the 9th column
#'   \item All features contain an ID or Name field
#' }
#' 
#' \preformatted{
#' Chr1   .   mRNA   3631   5899   .   +   .   ID=AT1G01010.1
#' Chr1   .    CDS   3760   3913   .   +   .   ID=AT1G01010.1.CDS-1;Parent=AT1G01010.1
#' Chr1   .    CDS   3996   4276   .   +   .   ID=AT1G01010.1.CDS-2;Parent=AT1G01010.1
#' }
#'
#' Expected extension: *.gff
#'
#' @section Genome input:
#'
#' This must be a fasta file (extension 'fna', for Fasta Nucleic Acid). The
#' header must contain sequence ids that match those of the GFF.
#'
#' @section Input synteny maps:
#'
#' Absolute path to a directory containing one synteny map for each species that
#' will be compared. each synteny map should consist of a single file named
#' according to the pattern "<query>.vs.<target>.syn", for example,
#' "arabidopsis_thaliana.vs.arabidopsis_lyrata.tab". these files should contain
#' the following columns:
#'
#' \enumerate{
#'   \item query contig name (e.g. chromosome or scaffold)
#'   \item query start position
#'   \item query stop position
#'   \item target contig name
#'   \item target start position
#'   \item target stop position
#'   \item score (not necessarily used)
#'   \item strand relative to query
#' }
#'
#' Example:
#'
#' \preformatted{
#' chr2   193631   201899   tchr2   193631   201899   100   +
#' chr2   225899   235899   tchr2   201999   202999   100   +
#' chr1   5999     6099     tchr1   6099     6199     100   +
#' chr1   5999     6099     tchr1   8099     8199     100   +
#' chr1   17714    18714    tchr2   17714    18714    100   +
#' chr2   325899   335899   tchr2   301999   302999   100   +
#' }
#'
#' a synteny map like this can be created using a whole genome synteny program,
#' such as satsuma (highly recommended). building a single synteny map requires
#' hundreds of cpu hours, so it is best done on a cluster. an example pbs script
#' is provided, see src/satsuma.pbs.
#'
#' Expected filename format: <query_sciname>.vs.<target_sciname>.syn
#'
#'
#' @section Input Tree:
#'
#' Absolute path to a newick format file specifying the topology of the species
#' tree. It must contain all species used in the pipeline AND NO OTHERS (I may
#' relax this restriction later).
#'
#' NOTE: There must be no spaces in the species names.
#'
#' Here is an example tree:
#'
#' (Brassica_rapa,(Capsella_rubella,(Arabidopsis_lyrata,Arabidopsis_thaliana)));
#'
#' @section Synder parameters:
#'
#' See documentation in \code{synder}
#'
#' @section Fagin parameters:
#'
#' \describe{
#' 
#'   \item{PROT2PROT_PVAL}{
#'     default=0.05 - Base p-value cutoffs. These will be ladjusted for multiple
#'     testing query protein versus target gene alignments.
#'   }
#' 
#'   \item{PROT2ALLORF_PVAL}{
#'     default=0.05 - query protein versus all SI translated ORFs.
#'   }
#' 
#'   \item{PROT2TRANSORF_PVAL}{
#'     default=0.05 - query protein versus translated ORFs from spliced transcripts
#'   }
#' 
#'   \item{DNA2DNA_PVAL}{
#'     default=0.05 - query genes versus entire SI (nucleotide match)
#'   }
#' 
#'   \item{PROT2PROT_NSIMS}{default=1000 - number of simulations}
#'   \item{PROT2ALLORF_NSIMS}{default=1000 - number of simulations}
#'   \item{PROT2TRANSORF_NSIMS}{default=1000 - number of simulations}
#' 
#'   \item{DNA2DNA_MAXSPACE}{default=1e8 - Maximum value of m*n that will be searched}
#' 
#'   \item{INDEL_THRESHOLD}{default=0.25 - Ratio of search interval to query
#'     interval below which an indel is called}
#' 
#' }
#' 
#' @docType package
#' @name fagin
NULL



#' Run a Fagin analysis
#'
#' @export
#' @param con A fagin_config object
#' @return An rmonad object containing all results
run_fagin <- function(con){
  {

    "Set random seed for the analysis, the choice of 210 is arbitrary. The
    random seed mainly affects the p-value estimates for alignments."

    set.seed(210)

  } %__% {

    "Record the date, system info, and installed packages. Of particular
    importance are the versions of `synder` and `rmonad`."

    devtools::session_info()
  
  } %v__% {

    "Create the archival directory"

    dir.create(con@archive)

  } %__%
  primary_data(con)      %>_% archive_1(con@archive) %>>%
  secondary_data(con)    %>_% archive_2(con@archive) %>>%
  tertiary_data(con)     %>_% archive_3(con@archive) %>>%
  determine_labels(con)  %>_% archive_4(con@archive) %>>%
  determine_origins(con) %>_% archive_5(con@archive) %>_%
                              archive_rmonad(con@archive)
}
