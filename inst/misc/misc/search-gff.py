#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict

def parser():
    parser = argparse.ArgumentParser(
        description=\
        '''
            Given a gff file and a synteny file, print the overlapping
            portions. The algorithm used is extremely slow (linear lookup), so
            is recommended only for GFFs with a few entries. Basically, it is
            useful for diagnostics.
        '''
    )
    parser.add_argument(
        'gff',
        help="GFF file to search",
        type=argparse.FileType('r')
    )
    parser.add_argument(
        'syn',
        help="Synteny file",
        type=argparse.FileType('r')
    )
    parser.add_argument(
        '-l', '--left-flank',
        help="Length of left flank to include in overlaps",
        type=int,
        default=0
    )
    parser.add_argument(
        '-r', '--right-flank',
        help="Length of right flank to include in overlaps",
        type=int,
        default=0
    )
    parser.add_argument(
        '-p', '--gff-pattern',
        metavar='PAT',
        help="Use the GFF row only if it matches PAT (a regular expression)"
    )
    parser.add_argument(
        '-f', '--feature',
        metavar='FEAT',
        help="Only search gff entries with FEAT as the 3rd field (e.g. mRNA, CDS, gene)"
    )
    args = parser.parse_args()

    return(args)

if __name__ == '__main__':
    args = parser()
    syn = defaultdict(list)
    for line in args.syn.readlines():
        terms = line.strip().split('\t')
        chrid = terms[0]
        a = int(terms[1])
        b = int(terms[2])
        syn[chrid].append((a, b, line.strip()))
    for line in args.gff.readlines():
        if args.gff_pattern and not re.search(args.gff_pattern, line):
            continue
        terms = line.strip().split('\t')
        chrid = terms[0]
        feat  = terms[2]
        start = int(terms[3])
        stop  = int(terms[4])
        seqid = terms[8]
        g1 = max(0, start - args.left_flank)
        g2 = stop + args.right_flank
        if args.feature and not args.feature == feat:
            continue
        entries = []
        for s1, s2, synline in syn[chrid]:
            if g1 <= s2 and g2 >= s1:
                rel_lo = s1 - start
                rel_hi = s2 - start
                if rel_hi > 0 and rel_lo < (stop - start):
                    overlaps = '$'
                else:
                    overlaps = ''
                out = "(%s,%s)%s " % (rel_lo, rel_hi, overlaps) + synline
                entries.append((rel_lo, out))
        print('>>>|%s|%s|%s|%s|%s' % (chrid, feat, str(g1), str(g2), seqid))
        for pos, line in sorted(entries):
            print(line)
