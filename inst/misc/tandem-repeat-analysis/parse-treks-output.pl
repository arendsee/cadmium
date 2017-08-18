#!/usr/bin/env perl
# repeat not found in sequence >AT1G02320.1
# repeat not found in sequence >AT1G02350.1
# repeat not found in sequence >AT1G02490.1
# repeat not found in sequence >AT1G03200.1
# repeat not found in sequence >AT1G14642.1
# >AT1G15840.1
# Length: 6 residues - nb: 10  from  41 to 100 - Psim:0.72 region Length:60 
# GEG-G-G---
# GEG-T-S---
# GEG-G-G-G-
# G-G-D-G-TK
# GGG-D-GIS-
# G-G-GHG-D-
# GLGCS-G-G-
# G-G-D-G-TK
# G-G-G-RRG-
# -DGLG-R-G-
# **********************

use strict;
use warnings;
use v5.10;

use Getopt::ArgParse;

my $parser = Getopt::ArgParse->new_parser(
    prog        => 'parse-treks-output.pl',
    description => 'Parse MSA and tab-delimited data from T-Reks '
);
$parser->add_arg('input',  required => 1);
$parser->add_arg('tabout', required => 1);
$parser->add_arg('msaout', required => 1);

my $args = $parser->parse_args();

open(my $INPUT, "<", $args->input)
    or die "Can't open < $args->input: $!";

open(my $TABOUT, ">", $args->tabout)
    or die "Can't open > $args->tabout: $!";

open(my $MSAOUT, ">", $args->msaout)
    or die "Can't open > $args->msaout: $!";

my $seqid;
my $msa_count;
my $in_msa = 0;
while(<$INPUT>){

    # Find proteins with no repeats found
    if(/^repeat not found in sequence >(\S+)$/) {
        print $TABOUT "$1_1\t0\t0\t0\t0\t0\t0\n";
    }

    # Parse out the sequence name of proteins with repeats
    elsif(/^>(\S+)$/){
        # Number of MSA encountered so far under current seqid
        $msa_count = 0; # reset the MSA count 
        $seqid = $1;    # storing the seqid for future use
    }

    # Parser for extracting line of info for identified motif
    # Length: 6 residues - nb: 10  from  41 to 100 - Psim:0.72 region Length:60 
    elsif(/^\s*Length:\s+(\d+)\s+residues\s+-\s+nb:\s+(\d+)\s+from\s+(\d+)\s+to\s+(\d+)\s+-\s+Psim:(..\d\d)\d*\s+region\s+Length:(\d+)\s*$/){
        $msa_count++;
        print $TABOUT "${seqid}_$msa_count\t$1\t$2\t$3\t$4\t$5\t$6\n";
        print $MSAOUT ">${seqid}_$msa_count\n";
        $in_msa = 1
    }

    # Match the end of an alignment
    elsif(/\*\*\*/){
        $in_msa = 0;
    }

    # Parser for writing to MSA file
    elsif($in_msa){
        print $MSAOUT $_;
    }

}

close $INPUT;
close $TABOUT;
close $MSAOUT;

exit;
