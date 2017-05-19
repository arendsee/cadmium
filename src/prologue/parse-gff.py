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

  Several constraints on the output entries may be enacted. These include the
  options --hasParent, --required, --uniqueTag, --requiredTag. All of these are
  tested on the output entries, after the --select and --reduce steps.

  If --uniqueTag=TAGS is set, any tag in TAGS is required to have unique values
  for every entry in the GFF. If --forceUnique is set as well, the parser
  tries to make non-unique values unique. It does this be appending
  "<TYPE>_<START>" to all values in the type that was not unique. For example,
  if ID is required to be unique, but a GFF file associates the same ID with
  all CDS in a given mRNA, then ID entries such as "foo123" will be replaced
  with "foo123_CDS_1066" (where 1066 is the starting position). If this
  addition does not result in unique ids, an error is still emitted. 

USAGE:
  parse-gff.py [--select=SELECT] [--reduce=REDUCE]
               [--required=TYPES] [--hasParent=KIDS] 
               [--requiredTag=TAGS] [--uniqueTag=TAGS] [--forceUnique=TAGS] 
               [--split] [--mapid] [--swapid] [--strict]
               [--untagged-ignore] [--untagged-name=TAG] GFFFILE
  parse-gff.py (-h | --help)
  parse-gff.py (--version)

OPTIONS:
  -h, --help           show this help message and exit
  --select=SELECT      filter 3rd column with these comma-delimited patterns 
  --reduce=REDUCE      reduce the 9th column to these comma-delimited tags
  --hasParent=KIDS     list of types that are required to have Parent tags
  --required=TYPES     list of types that MUST be present
  --uniqueTag=TAGS     list of tags that MUST be unique (but may be missing)
  --forceUnique=TAGS   try to fix non-unique tags in --uniqueTag (see note)
  --requiredTag=TAGS   list of tags that must be present
  --split              split attribute column by tag (in REDUCE order)
  --mapid              map parent ID to parent Name
  --swapid             if no "Name" field is present, use ID instead
  --strict             fail if anything looks fishy (see notes above)
  --untagged-ignore    ignore attributes with no tag
  --untagged-name=TAG  give attributes with no tag the name TAG [default=None]

EXAMPLES:
1) Extract elements from a GFF, require all types are represented:
     parse-gff.py --select=mRNA,exon --strict foo.gff

2) Extract CDS from a GFF, require they all have a parent listed, that they
   have a unique ID, then include only the ID and Parent tags
     parse-gff.py --select=CDS --strict --hasParent=CDS --reduce=ID,Parent \\
                  --uniqTags=ID foceUnique=ID foo.gff
"""

import sys
import signal
from collections import defaultdict
from docopt import docopt

if sys.version_info[0] != 3 or sys.version_info[1] < 5:
    sys.exit("This script requires Python 3.5+")


# Global variable, the name of the input file
GLOBAL_GFF = None
# It is used only here, to generate error messages
def err(msg):
    header = "%s in %s" % (Colors.error("ERROR:"), GLOBAL_GFF)
    print(header, file=sys.stderr, end="\n  ")
    print(msg, file=sys.stderr)
    sys.exit(1)


class Error:

    def handle_no_output():
        err("No output resulting from parsing '%s'")

    def handle_untagged_attribute(attributes):
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

    def handle_invalid_gff(row):
        msg = "Bad GFF, must be TAB-delimited with 9 columns\n" 
        msg += "  Offending line (TABs substituted for '|'):\n"
        msg += "  %s" % '|'.join(row)
        err(msg)

    def handle_missing_id():
        err("Bad GFF, 9th column must have ID tag")


    def handle_missing_selected_item(expstr, missing, filename):
        msg = "expected the types {} in '{}', but {} missing"
        msg = msg.format(expstr, filename, missing);
        err(msg)

    def handle_duplicate_attribute(key, row):
        msg = "Tag '{}' is duplicated in attribute column:\n  > {}".format(key,row)
        err(msg)

    def handle_missing_parent(entry):
        msg  = "Entries or type '{}' are required to have a Parent=<ID> entry\n"
        msg += "in the 9th column. Offending line (TAB replace with '|'):\n  > {}"
        msg = msg.format(entry.row[2], entry.toString(sep="|"))
        err(msg)


    def handle_missing_required_type(missing, all_seen):
        msg  = "Missing required entry type %s, actual types in file are:\n" % str(missing)
        msg += '\n'.join(["  * %s" % x for x in all_seen])
        err(msg)

    def handle_unique_tag_violation(tag):
        msg = "The tag '%s' is required to be unique, but isn't" % tag
        err(msg)

    def handle_required_tags_violation(tag, e):
        msg  = "The required tag '%s' is missing\n"
        msg += "first offending line (TABs replaced with '|'):\n%s\n"
        msg = msg % (tag, e.toString(sep="|"))
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
                    Error.handle_untagged_attribute(row[8].split(";"))
            if(not keepers or t in keepers):
                if t in self.attr:
                    Error.handle_duplicate_attribute(t, row[8])
                else:
                    self.attr[t] = v

    def toString(self, tags=None, split=False, use_ids=False, sep="\t"):
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

        out = sep.join(self.row[0:8] + [nine])

        return out

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
            Error.handle_invalid_gff(row)

        entries.append(Entry(row, keepers, **kwargs))

    return entries


def mapid(entries):
    idmap = dict()
    for entry in entries:
        try:
            ID = entry.attr['ID']
        except KeyError:
            Error.handle_missing_id()
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


def tag_is_unique(tag, entries):
    ''' Test whether each value associated with given tag is unique '''
    values = [e.attr[tag] for e in entries if tag in e.attr]
    is_unique = len(values) == len(set(values))
    return is_unique

def _force_unique_tag_type(tag, typ, entries):
    for e in entries:
        if e.row[2] == typ:
            e.attr[tag] = "%s_%s-%s" % (e.attr[tag], e.row[2], e.row[3])
    return entries

def force_unique(tag, entries):
    dic = defaultdict(list)
    for e in entries:
        dic[e.row[2]].append(e)
    for k,v in dic.items():
        if not tag_is_unique(tag, v):
            entries = _force_unique_tag_type(tag, k, entries)
    return entries


def condition_unique_tags(tags, entries):
    tags = set(tags)
    for tag in tags:
        if not tag_is_unique(tag, entries):
            Error.handle_unique_tag_violation(tag)


def condition_all_selected_in_file(selected, entries, filename):
    obs_types = set([e.row[2] for e in entries])
    exp_types = set([s for s in selected])
    if exp_types and obs_types != exp_types:
        expstr = "[%s]" % (','.join(exp_types))
        missing = "[%s]" % (','.join(exp_types - obs_types))
        Error.handle_missing_selected_item(expstr, missing, filename)


def condition_children_have_parents(kids, entries):
    kids = set(kids)
    for e in entries:
        if (e.row[2] in kids) and not ("Parent" in e.attr):
            Error.handle_missing_parent(e)


def condition_required_types(required, entries):
    types = set([e.row[2] for e in entries])
    required = set(required)
    if (required - types):
        Error.handle_missing_required_type(required - types, types)


def condition_required_tags(tags, entries):
    tags = set(tags)
    for tag in tags:
        for e in entries:
            if not tag in e.attr:
                Error.handle_required_tags_violation(tag, e)


if __name__ == '__main__':

    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    args = docopt(__doc__, version='0.2.0')

    GLOBAL_GFF = args['GFFFILE']

    for tag in ["--reduce", "--select", "--hasParent", "--required", "--uniqueTag", "--requiredTag"]:
        if args[tag]:
            args[tag] = args[tag].split(",")

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

        if args["--forceUnique"] and args["--uniqueTag"]:
            for tag in args["--forceUnique"]:
                entries = force_unique(tag, entries)


        # ---------------------------- CONDITIONS ----------------------------
        if args["--hasParent"]:
            condition_children_have_parents(args["--hasParent"], entries)

        if args["--required"]:
            condition_required_types(args["--required"], entries)

        #  if args["--uniqueTag"]:
        #      condition_unique_tags(args["--uniqueTag"], entries)

        if args["--requiredTag"]:
            condition_required_tags(args["--requiredTag"], entries)

        if args["--strict"] and args["--select"]:
            condition_all_selected_in_file(
                args["--select"],
                entries,
                args["GFFFILE"]
            )

        if args["--strict"] and not entries:
            Error.handle_no_output(args["GFFFILE"])
        # --------------------------------------------------------------------

        for entry in entries:
            s = entry.toString(args["--reduce"], attr_join, args["--swapid"])
            print(s)

    sys.exit(0)
