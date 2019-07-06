#!/usr/bin/env python3

from argparse import ArgumentParser


GAP_MAIN_LOOP = \
"""for d in [{},{}..{}] do
{{}}
od;"""

GAP_PRIMITIVE_GROUP = \
"""  for n in [1..NrPrimitiveGroups(d)] do
    g := PrimitiveGroup(d, n);
    gens := GeneratorsOfGroup(g);
    Print(gens, "\\n");
  od;"""

GAP_SPECIAL_GROUP = \
"""  g := {}(d);
  gens := GeneratorsOfGroup(g);
  Print(gens, "\\n");"""


if __name__ == '__main__':
    parser = ArgumentParser(
        description="Obtain generators of known permutation groups")

    parser.add_argument('group_class', metavar='GROUP_CLASS',
      choices=['symmetric', 'cyclic', 'alternating', 'dihedral', 'primitive'],
      help="class of groups for which to fetch generators")

    parser.add_argument('--min-degree', default=1,
        help="minimum group degree")

    parser.add_argument('--max-degree', required=True,
        help="maximum group degree")

    parser.add_argument('--degree-step', default=1,
        help="group degree step")

    ns = parser.parse_args()

    gap_cmd = GAP_MAIN_LOOP.format(
        ns.min_degree, ns.min_degree + ns.degree_step, ns.max_degree)

    if ns.group_class == 'primitive':
        gap_cmd = gap_cmd.format(GAP_PRIMITIVE_GROUP)
    else:
        gap_extract_gens = GAP_SPECIAL_GROUP.format({
            'symmetric':   'SymmetricGroup',
            'cyclic':      'CyclicGroup',
            'alternating': 'AlternatingGroup',
            'dihedral':    'DihedralGroup'
        }[ns.group_class])

        gap_cmd = gap_cmd.format(gap_extract_gens)

    print(gap_cmd)
