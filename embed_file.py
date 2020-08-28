#!/usr/bin/env python3

from argparse import ArgumentParser


def to_array(data):
    return '{' + ','.join(hex(ord(c)) for c in data) + '}'


if __name__ == '__main__':
    parser = ArgumentParser(description="Convert text file into embeddable C source file")

    parser.add_argument('resource', metavar='RESOURCE')
    parser.add_argument('outfile', metavar='outfile')
    parser.add_argument('symbol', metavar='SYMBOL')

    ns = parser.parse_args()

    res = ns.resource
    out = ns.outfile
    sym = ns.symbol

    with open(res, 'r') as f:
        data = f.read().rstrip()

    with open(out, 'w') as f:
        f.write('#include <stddef.h>\n')
        f.write(f'char const {sym}[] = {to_array(data)};\n')
        f.write(f'size_t const {sym}_len = sizeof({sym});\n')
