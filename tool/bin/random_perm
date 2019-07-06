#!/usr/bin/env python3

import sys

from argparse import ArgumentParser
from random import shuffle
from typing import List


def rand_perm(n: int) -> List[int]:
    """
    Generate a random permutation.

    :param n: determines the permutation domain `{1,...,n}`
    :returns: a partial permutation represented by a list of index mappings
    """
    perm = [x for x in range(1, n + 1)]
    shuffle(perm)

    return perm


if __name__ == '__main__':
    parser = ArgumentParser(description="Generate random permutations.")

    parser.add_argument('N',
        help="maximum domain element n of domain {1,...,n}")

    parser.add_argument('-n', default=1, type=int,
        help="number of partial permutations to be generated")

    parser.add_argument('--sort', action='store_true',
        help="sort result alphabetically")

    ns = parser.parse_args()

    res = []
    for _ in range(ns.n):
        perm = rand_perm(int(ns.N))
        strlist = ', '.join(map(str, perm))
        res.append("Perm({{{}}})".format(strlist))

    if ns.sort:
        res = sorted(res)

    print(',\n'.join(res))
