#!/usr/bin/env python3

import os
import sys

from argparse import ArgumentParser, RawTextHelpFormatter
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile


USAGE = \
"""get_perm_groups.py [-h]
                          [--min-degree MIN_DEGREE]
                           --max-degree MAX_DEGREE
                          [--degree-step DEGREE_STEP]
                          [--order-limit ORDER_LIMIT]
                          CLASS"""

DESC = "Obtain generators of known permutation groups"

GAP_CMD = \
"""for d in [{},{}..{}] do
     {{}}
   od;"""

GAP_ORDER_SKIP = "if Order(g) > {} then continue; fi;"

GAP_GET_GENERATORS = "gens := GeneratorsOfGroup(g);"

GAP_PRINT = 'Print("degree:", d, ",order:", Order(g), ",gens:", gens, "\\n");'

GAP_PRIMITIVE_GROUP = \
"""for n in [1..NrPrimitiveGroups(d)] do
    g := PrimitiveGroup(d, n);
    {}
    {}
    {}
  od;"""

GAP_SPECIAL_GROUP = \
"""g := {}(d);
   {{}}
   {{}}
   {{}}"""


def gap_script(group_class, min_degree, max_degree, degree_step, order_limit):
    gap_cmd = GAP_CMD.format(min_degree, min_degree + degree_step, max_degree)

    if order_limit is not None:
        gap_order_skip = GAP_ORDER_SKIP.format(order_limit)
    else:
        gap_order_skip = ""

    if group_class == 'primitive':
        gap_group = GAP_PRIMITIVE_GROUP
    else:
        gap_group = GAP_SPECIAL_GROUP.format({
            'symmetric':   'SymmetricGroup',
            'cyclic':      'CyclicGroup',
            'alternating': 'AlternatingGroup',
            'dihedral':    'DihedralGroup'
        }[group_class])

    gap_group = gap_group.format(gap_order_skip, GAP_GET_GENERATORS, GAP_PRINT)

    return gap_cmd.format(gap_group)


if __name__ == '__main__':
    def formatter_class(prog):
        return RawTextHelpFormatter(prog, max_help_position=80)

    parser = ArgumentParser(
        usage=USAGE,
        description=DESC,
        formatter_class=formatter_class)

    parser.add_argument('group_class', metavar='CLASS',
      choices=['symmetric', 'cyclic', 'alternating', 'dihedral', 'primitive'],
      help="groups class")

    parser.add_argument('--min-degree', default=1,
        help="minimum group degree, default is %(default)d")

    parser.add_argument('--max-degree', required=True,
        help="maximum group degree")

    parser.add_argument('--degree-step', default=1,
        help="group degree step, default is %(default)d")

    parser.add_argument('--order-limit', default=2**63-1,
        help=str("group order limit, default is %(default)d"))

    ns = parser.parse_args()

    script = gap_script(ns.group_class,
                        ns.min_degree,
                        ns.max_degree,
                        ns.degree_step,
                        ns.order_limit)

    with NamedTemporaryFile(mode='w+', delete=False) as f:
        f.write(script)

    p = Popen(['gap', '--nointeract', '-q', f.name], stdout=PIPE)
    out, _ = p.communicate()
    out = out.decode('ascii')[:-1]

    out = out.replace(' ', '')

    out_list = [c for c in out]

    i = 0
    while i < len(out_list):
        if out_list[i] == '\n' and out_list[i - 1] in ':),':
            out_list[i] = ''
        i += 1

    print(''.join(out_list))

    os.remove(f.name)
