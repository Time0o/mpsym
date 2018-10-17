#!/usr/bin/env python3

import argparse
import os
import re
import sys
from typing import List, Optional

from argparse import ArgumentParser

USAGE = ""

DESC = """\
Reconstructs valid GAP partial permutation construction statements from their
string representations passed to this program via stdin."""

PPERM = r'^((?:\[\d+,\d+(?:,\d+)*\])*)((?:\(\d+(?:,\d+)*\))*)$'
PPERM_IDENTITY = r'^<identity partial perm on \[ (\d+(?:, \d+)*) \]>$'


def pperm_identity_constructor(domain: List[int]) -> str:
    return 'PartialPermOp((), {})'.format(domain)


def pperm_constructor(chains: Optional[List[List[int]]],
                      cycles: Optional[List[List[int]]]) -> str:
    return ''  # TODO


if __name__ == '__main__':
    # program name
    progname = os.path.basename(__file__)

    # parse arguments
    def formatter_class(prog):
        return argparse.RawTextHelpFormatter(prog, max_help_position=80)

    parser = ArgumentParser(description=DESC, formatter_class=formatter_class)

    parser.add_argument('-f', '--infile',
                        help="read from IN_FILE instead of stdin")

    help_in_place = "if used together with --in-file, modify input file in place"
    parser.add_argument('-i', '--in-place',
                        action='store_true',
                        help=help_in_place)

    parser.add_argument('-o', '--outfile',
                        help="write to OUT_FILE instead of stdout")

    parser.add_argument('-s', '--sort', action='store_true',
                        help="sort result alphabetically")

    args = parser.parse_args()

    if args.in_place:
        if args.infile is None:
            parser.print_usage(sys.stderr)

            fmt = "{}: error: --in-place may only be used in combination with --in-file"
            print(fmt.format(progname), file=sys.stderr)

            sys.exit(1)

        if args.outfile is not None:
            parser.print_usage(sys.stderr)

            fmt = "{}: error: --in-place may not be used in combination with --out-file"
            print(fmt.format(progname), file=sys.stderr)

            sys.exit(1)

    # parse input
    matcher_pperm_identity = re.compile(PPERM_IDENTITY)
    matcher_pperm = re.compile(PPERM)

    result = []

    istream = open(args.infile, 'r') if args.infile is not None else sys.stdin

    for line in istream:
        line = line.rstrip()

        match_pperm_identity = matcher_pperm_identity.match(line)
        if match_pperm_identity is None:

            match_pperm = matcher_pperm.match(line)
            chains, cycles = match_pperm.groups()

            if match_pperm is None or not (chains or cycles):
                fmt = "{}: error: failed to parse input line:\n{}"
                print(fmt.format(progname, line), file=sys.stderr)

                sys.exit(1)

            chains = [chain.split(',') for chain in chains[1:-1].split('][')]
            cycles = [cycles.split(',') for cycles in cycles[1:-1].split(')(')]

            result.append(pperm_constructor(chains, cycles))

        else:
            domain = match_pperm_identity.group(1).split(',')

            result.append(pperm_identity_constructor(domain))

    if istream != sys.stdin:
        istream.close()

    if args.sort:
        result = sorted(result)

    result = ',\n'.join(result)

    if args.in_place:
        with open(args.infile, 'w', encoding='UTF-8') as outfile:
            outfile.write(result + '\n')

    elif args.outfile is not None:
        with open(args.outfile, 'w', encoding='UTF-8') as outfile:
            outfile.write(result + '\n')

    else:
        print(result)
