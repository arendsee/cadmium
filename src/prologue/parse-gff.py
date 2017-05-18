#!/usr/bin/env python3

"""
merge-gff.py - prepare GFF file Fagin analysis

DESCRIPTION:
  Fagin requires detailed information about gene models. It requires clear
  links between the parent and child entries: a "gene" entry is the parent of
  multiple "mRNA"; an "mRNA" entry is the parent of many "exon", "CDS", and UTR
  entries. Every entry also needs a unique ID. However, GFF files are not well
  standardized.

  `merge-gff.py` will attempt to uncover all parent-child relations and set
  unique ids, if it fails, it will report what went wrong and what you need to
  do to fix it.

  The `--strict` option will enforce the following conditions:
    1. All values specied in `--select` must appear in file
    2. At least one row is output

USAGE:
  parse-gff.py [--select=SELECT] [--reduce=REDUCE] [--split] [--mapid]
               [--swapid] [--strict] [--untagged-ignore]
               [--untagged-name=TAG] GFFFILE
  parse-gff.py (-h | --help)
  parse-gff.py (--version)

OPTIONS:
  -h, --help          show this help message and exit
  --select=SELECT     filter 3rd column with these comma-delimited patterns 
  --reduce=REDUCE     reduce the 9th column to these comma-delimited tags
  --split             split attribute column by tag (in REDUCE order)
  --mapid             map parent ID to parent Name
  --swapid            if no "Name" field is present, use ID instead
  --strict            fail if anything looks fishy (see notes above)
  --untagged-ignore   ignore attributes with no tag
  --untagged-name=TAG give attributes with no tag the name TAG [default=None]

EXAMPLES:
1) Extract elements from a GFF, require all types are represented:
   parse-gff.py --select=mRNA,exon --strict foo.gff
"""

import sys
import signal
from docopt import docopt

if sys.version_info[0] != 3 or sys.version_info[1] < 5:
    sys.exit("This script requires Python 3.5+")


class Error:

    def handle_untagged_attribute_error(attributes):
        elements = []
        for s in attributes:
            if len(s.split('=')) == 2:
                elements.append(s) 
            else:
                elements.append(Colors.error(s))
        msg = \
"""inappropriately formatted attribute column
  expected /<tag>=<value>(;<tag>=<value>)*/
  offending attribute entry:
    > """ + ";".join(elements) + Colors.emphasis("\nPossible solutions:") + """
  1) Remove untagged attributes (--untagged-ignore)
  2) If there is only one unnamed attribute, give it the name TAG
     (--untagged-name=TAG) If there are multiple unnamed attributes,
     this will still die
"""
        err(msg)

class Colors:
    OFF          = chr(27) + '[0;0m'
    RED          = chr(27) + '[0;31m'
    GREEN        = chr(27) + '[0;32m'
    YELLOW       = chr(27) + '[0;33m'
    MAGENTA      = chr(27) + '[0;35m'
    CYAN         = chr(27) + '[0;36m'
    WHITE        = chr(27) + '[0;37m'
    BLUE         = chr(27) + '[0;34m'
    BOLD_RED     = chr(27) + '[1;31m'
    BOLD_GREEN   = chr(27) + '[1;32m'
    BOLD_YELLOW  = chr(27) + '[1;33m'
    BOLD_MAGENTA = chr(27) + '[1;35m'
    BOLD_CYAN    = chr(27) + '[1;36m'
    BOLD_WHITE   = chr(27) + '[1;37m'
    BOLD_BLUE    = chr(27) + '[1;34m'

    def _color(s, color):
        if sys.stdout.isatty():
            s = color + s + Colors.OFF
        return s

    def error(s):
        return Colors._color(s, Colors.BOLD_RED)

    def good(s):
        return Colors._color(s, Colors.BOLD_GREEN)

    def warn(s):
        return Colors._color(s, Colors.BOLD_YELLOW)

    def emphasis(s):
        return Colors._color(s, Colors.BOLD_WHITE)


def err(msg):
    print(Colors.error("ERROR:"), file=sys.stderr, end=" ")
    print(msg, file=sys.stderr)
    sys.exit(1)


class Entry:
    def __init__(self, row, keepers, untagged_ignore=False, untagged_name=None):
        self.row = row[0:8]
        self.attr = dict()
        for attribute in row[8].split(';'):
            try:
                t,v = attribute.split('=')
            except ValueError:
                if untagged_name:
                    t = untagged_name
                    v = attribute
                elif untagged_ignore:
                    continue 
                else:
                    Error.handle_untagged_attribute_error(row[8].split(";"))
            if(not keepers or t in keepers):
                self.attr[t] = v

    def print_(self, tags=None, split=False, use_ids=False):
        if(tags):
            if use_ids and 'Name' in tags and not 'Name' in self.attr:
                self.attr['Name'] = self.get_attr('ID')
            if split:
                nine = '\t'.join([self.get_attr(k) for k in tags])
            else:
                if len(tags) == 1:
                    nine = self.get_attr(tags[0])
                else:
                    nine = ';'.join(['%s=%s' % (k, self.get_attr(k)) for k in tags])
        else:
            nine = ';'.join(['%s=%s' % (k,v) for k,v in self.attr.items()])

        out = '\t'.join(self.row[0:8] + [nine])

        print(out)

    def get_attr(self, tag):
        try:
            return self.attr[tag]
        except KeyError:
            return '-'


def rowgen(gff, keepers, **kwargs):
    entries = []
    for line in gff.readlines():
        if(line[0] == '#'):
            continue

        row = line.rstrip().split('\t')

        if(len(row) != 9):
            msg = "Bad GFF, must be TAB-delimited with 9 columns\n" 
            msg += "  Offending line (TABs substituted for '|'):\n"
            msg += "  %s" % '|'.join(row)
            err(msg)

        entries.append(Entry(row, keepers, **kwargs))

    return entries


def mapid(entries):
    idmap = dict()
    for entry in entries:
        try:
            ID = entry.attr['ID']
        except KeyError:
            err("Bad GFF, 9th column must have ID tag")

        try:
            Name = entry.attr['Name']
        except KeyError:
            Name = entry.attr['ID']

        idmap[ID] = Name
    return(idmap)


def parent_id2name(entries, idmap):
    for entry in entries:
        try:
            entry.attr['Parent'] = idmap[entry.attr['Parent']]
        except KeyError:
            pass
        except TypeError:
            pass


def condition_all_selected_in_file(selected, entries, filename):
    obs_types = set([e.row[2] for e in entries])
    exp_types = set([s for s in selected])
    if exp_types and obs_types != exp_types:
        msg = "expected the types {} in '{}', but {} missing"
        expstr = "[%s]" % (','.join(exp_types))
        missing = "[%s]" % (','.join(exp_types - obs_types))
        msg = msg.format(expstr, filename, missing)
        err(msg)


if __name__ == '__main__':

    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    args = docopt(__doc__, version='0.2.0')

    if args["--reduce"]:
        args["--reduce"] = args["--reduce"].split(",")

    if args["--select"]:
        args["--select"] = args["--select"].split(",")

    keepers = set()
    if(args["--reduce"]):
        keepers = set(args["--reduce"])
        if(args["--mapid"] or args["--swapid"]):
            keepers.update(['ID', 'Name'])

    with open(args["GFFFILE"], "r") as gff:

        entries = rowgen(
            gff,
            keepers,
            untagged_ignore=args["--untagged-ignore"],
            untagged_name=args["--untagged-name"]
        )

        if args["--mapid"]:
            idmap = mapid(entries)
            parent_id2name(entries, idmap)

        attr_join = '\t' if args["--split"] else None

        if args["--select"]:
            entries = [e for e in entries if e.row[2] in args["--select"]]


        if args["--strict"]:
            condition_all_selected_in_file(
                args["--select"],
                entries,
                args["GFFFILE"]
            )

            if not entries:
                err("No output resulting from parsing '%s'" % args["GFFFILE"])

        for entry in entries:
            entry.print_(args["--reduce"], attr_join, args["--swapid"])

    sys.exit(0)
