#!/usr/bin/env python3

import argparse
import sys

if sys.version_info[0] != 3 or sys.version_info[1] < 5:
    sys.exit("This script requires Python 3.5+")

def err(msg):
    sys.exit(msg)


class Entry:
    def __init__(self, row, keepers):
        self.row = row[0:8]
        self.attr = dict()
        for s in row[8].split(b';'):
            t,v = s.split(b'=')
            if(not keepers or t in keepers):
                self.attr[t] = v

    def print_(self, tags=None, split=False, use_ids=False):
        if(tags):
            if use_ids and b'Name' in tags and not b'Name' in self.attr:
                self.attr[b'Name'] = self.get_attr(b'ID')
            if split:
                nine = b'\t'.join([self.get_attr(k) for k in tags])
            else:
                if len(tags) == 1:
                    nine = self.get_attr(tags[0])
                else:
                    nine = b';'.join([b'%s=%s' % (k, self.get_attr(k)) for k in tags])
        else:
            nine = b';'.join([b'%s=%s' % (k,v) for k,v in self.attr.items()])

        out = b'\t'.join(self.row[0:8] + [nine])

        print(out.decode())

    def get_attr(self, tag):
        try:
            return self.attr[tag]
        except KeyError:
            return b'-'


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'gfffile',
        help="GFF input file",
        type=argparse.FileType('r')
    )
    parser.add_argument(
        '-s', '--select',
        help="select where 3rd column matches one of these (1+)",
        nargs='+'
    )
    parser.add_argument(
        '-r', '--reduce',
        help="reduce the attribute column (9th) to these tags",
        nargs='+'
    )
    parser.add_argument(
        '-p', '--split',
        help="split attribute column by tag (in REDUCE order)",
        action="store_true",
        default=False
    )
    parser.add_argument(
        '-m', '--mapids',
        help="map parent ID to parent Name",
        action="store_true",
        default=False
    )
    parser.add_argument(
        '-d', '--use-id-if-unnamed',
        help="If no Name field is present, use ID instead",
        action="store_true",
        default=False
    )
    args = parser.parse_args()
    return(args)


def rowgen(gff, keepers):
    entries = []
    for line in gff.readlines():
        line = line.encode()
        if(line[0] == ord('#')):
            continue

        row = line.rstrip().split(b'\t')

        if(len(row) != 9):
            err("Bad GFF, must be TAB-delimited with 9 columns")

        entries.append(Entry(row, keepers))

    return entries


def mapids(entries):
    idmap = dict()
    for entry in entries:
        try:
            ID = entry.attr[b'ID']
        except KeyError:
            err("Bad GFF, 9th column must have ID tag")

        try:
            Name = entry.attr[b'Name']
        except KeyError:
            Name = entry.attr[b'ID']

        idmap[ID] = Name
    return(idmap)


def parent_id2name(entries, idmap):
    for entry in entries:
        try:
            entry.attr[b'Parent'] = idmap[entry.attr[b'Parent']]
        except KeyError:
            pass
        except TypeError:
            pass


if __name__ == '__main__':
    args = parser()

    keepers = set()
    if(args.reduce):
        args.reduce = [s.encode() for s in args.reduce]
        keepers = set(args.reduce)
        if(args.mapids or args.use_id_if_unnamed):
            keepers.update([b'ID', b'Name'])

    entries = rowgen(args.gfffile, keepers)

    if(args.mapids):
        idmap = mapids(entries)
        parent_id2name(entries, idmap)

    attr_join = b'\t' if args.split else None

    if(args.select):
        selection = [s.encode() for s in args.select]
        entries = [e for e in entries if e.row[2] in selection]

    # Exit will succeed if there is something to print
    exit_status = 0 if bool(entries) else 1

    for entry in entries:
        entry.print_(args.reduce, attr_join, args.use_id_if_unnamed)

    sys.exit(exit_status)
