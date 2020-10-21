import pickle
import unittest
from copy import deepcopy
from itertools import cycle, permutations
from math import factorial
from random import sample
from textwrap import dedent

import pympsym as mp


class PermTest(unittest.TestCase):
    def setUp(self):
        self.dom = [0,1,2,3]

        self.p = mp.Perm([1,0,3,2])
        self.q = mp.Perm([1,2,0,3])

    def test_euqality(self):
        self.assertEqual(self.p, self.p)
        self.assertEqual(self.q, self.q)
        self.assertNotEqual(self.p, self.q)

    def test_getitem(self):
        self.assertEqual([self.p[x] for x in self.dom], [1,0,3,2])
        self.assertEqual([self.q[x] for x in self.dom], [1,2,0,3])

    def test_tosequence(self):
        self.assertEqual(list(self.p), [1,0,3,2])
        self.assertEqual(list(self.q), [1,2,0,3])

    def test_invert(self):
        self.assertEqual(~self.p, [1,0,3,2])
        self.assertEqual(~self.q, [2,0,1,3])

    def test_mul(self):
        self.assertEqual(self.p * self.q, [2,1,3,0])
        self.assertEqual(self.q * self.p, [0,3,1,2])
        self.assertEqual(list(self.p) * self.q, [2,1,3,0])
        self.assertEqual(self.q * list(self.p), [0,3,1,2])

    def test_bool(self):
        self.assertTrue(self.p)
        self.assertTrue(self.q)
        self.assertFalse(self.p * ~self.p)
        self.assertFalse(self.q * ~self.q)
        self.assertFalse(~self.p * self.p)
        self.assertFalse(~self.q * self.q)

    def test_degree(self):
        self.assertEqual(self.p.degree(), 4)
        self.assertEqual(self.q.degree(), 4)


class PermGroupTest(unittest.TestCase):
    def setUp(self):
        self.dom = [0,1,2,3]

        self.pg = mp.PermGroup([[1,0,3,2], [1,2,0,3]])

        self.pg_id = mp.PermGroup(4)

        self.pg_elems = {mp.Perm((0,1,2,3)),
                         mp.Perm((0,2,3,1)),
                         mp.Perm((0,3,1,2)),
                         mp.Perm((1,0,3,2)),
                         mp.Perm((1,2,0,3)),
                         mp.Perm((1,3,2,0)),
                         mp.Perm((2,0,1,3)),
                         mp.Perm((2,1,3,0)),
                         mp.Perm((2,3,0,1)),
                         mp.Perm((3,0,2,1)),
                         mp.Perm((3,1,0,2)),
                         mp.Perm((3,2,1,0))}

        self.pg_id_elems = {mp.Perm((0,1,2,3))}

    def test_len(self):
        self.assertEqual(len(self.pg), 12)
        self.assertEqual(len(self.pg_id), 1)

    def test_iter(self):
        self.assertEqual(set(self.pg), self.pg_elems)
        self.assertEqual(set(self.pg_id), self.pg_id_elems)

    def test_contains(self):
        for elem in permutations(self.dom):
          elem = mp.Perm(elem)

          if elem in self.pg_elems:
            self.assertTrue(elem in self.pg)
          else:
            self.assertFalse(elem in self.pg)

    def test_bool(self):
        self.assertTrue(self.pg)
        self.assertFalse(self.pg_id)

    def test_degree(self):
        self.assertEqual(self.pg.degree(), 4)
        self.assertEqual(self.pg_id.degree(), 4)

    def test_generators(self):
        self.assertEqual(self.pg.generators(), self.pg)

    def test_properties(self):
        self.assertFalse(self.pg.is_symmetric())
        self.assertTrue(self.pg.is_alternating())
        self.assertFalse(self.pg.is_symmetric())


class ArchGraphSystemTest(unittest.TestCase):
    HAEC_LUA = dedent(
        """
        local mpsym = require 'mpsym'

        local super_graph_clusters = mpsym.identical_clusters(4, 'SoC')
        local super_graph_channels = mpsym.linear_channels(super_graph_clusters, 'wireless')

        local proto_processors = mpsym.identical_processors(16, 'P')
        local proto_channels = mpsym.grid_channels(proto_processors, 'optical')

        return mpsym.ArchUniformSuperGraph:create{
          super_graph = mpsym.ArchGraph:create{
            directed = false,
            clusters = super_graph_clusters,
            channels = super_graph_channels
          },
          proto = mpsym.ArchGraph:create{
            directed = false,
            processors = proto_processors,
            channels = proto_channels
          }
        }
        """
    )

    def setUp(self):
        self.ag = mp.ArchGraphSystem.from_lua(self.HAEC_LUA)

        self.ag_orbit1 = [(0, 1, 2, 3),
                          (0, 4, 8, 12),
                          (3, 2, 1, 0),
                          (3, 7, 11, 15),
                          (12, 8, 4, 0),
                          (12, 13, 14, 15),
                          (15, 11, 7, 3),
                          (15, 14, 13, 12),
                          (48, 49, 50, 51),
                          (48, 52, 56, 60),
                          (51, 50, 49, 48),
                          (51, 55, 59, 63),
                          (60, 56, 52, 48),
                          (60, 61, 62, 63),
                          (63, 59, 55, 51),
                          (63, 62, 61, 60)]

        self.ag_orbit2 = [(0, 3, 12, 15),
                          (0, 12, 3, 15),
                          (3, 0, 15, 12),
                          (3, 15, 0, 12),
                          (12, 0, 15, 3),
                          (12, 15, 0, 3),
                          (15, 3, 12, 0),
                          (15, 12, 3, 0),
                          (48, 51, 60, 63),
                          (48, 60, 51, 63),
                          (51, 48, 63, 60),
                          (51, 63, 48, 60),
                          (60, 48, 63, 51),
                          (60, 63, 48, 51),
                          (63, 51, 60, 48),
                          (63, 60, 51, 48)]

    def test_num_processors(self):
        self.assertEqual(self.ag.num_processors(), 64)

    def test_num_channels(self):
        self.assertEqual(self.ag.num_channels(), 864)

    def test_automorphisms(self):
        autom = self.ag.automorphisms()
        self.assertEqual(len(autom), 8192)

    def test_representative(self):
        for orbit in [self.ag_orbit1, self.ag_orbit2]:
            for mapping in orbit:
                self.assertEqual(self.ag.representative(mapping), orbit[0])

    def test_orbit(self):
        for orbit in [self.ag_orbit1, self.ag_orbit2]:
            self.assertEqual(self.ag.orbit(orbit[0]), orbit)

    def test_from_nauty(self):
        vertices_super = 4
        adj_super = {0: [1], 1: [2], 2: [3]}

        vertices_proto = 16
        adj_proto = {0:  [1,4],  1:  [2,5],   2:  [3,6],   3:  [7],
                     4:  [5,8],  5:  [6,9],   6:  [7,10],  7:  [11],
                     8:  [9,12], 9:  [10,13], 10: [11,14], 11: [15],
                     12: [13],   13: [14],    14: [15]}

        ag_super = mp.ArchGraphSystem.from_nauty(vertices_super, adj_super, directed=False)
        ag_proto = mp.ArchGraphSystem.from_nauty(vertices_proto, adj_proto, directed=False)

        haec_nauty = mp.ArchUniformSuperGraph(ag_super, ag_proto)

        self.assertEqual(haec_nauty.automorphisms(), self.ag.automorphisms())

    def test_to_from_json(self):
        json = self.ag.to_json()
        ag_from_json = mp.ArchGraphSystem.from_json(json)
        self.assertEqual(ag_from_json.automorphisms(), self.ag.automorphisms())

    def test_expand_automorphisms(self):
        ag_expanded = self.ag.expand_automorphisms()
        self.assertEqual(ag_expanded.automorphisms(), self.ag.automorphisms())

    def test_deepcopy(self):
        ag_deepcopy = deepcopy(self.ag)
        self.assertEqual(ag_deepcopy.automorphisms(), self.ag.automorphisms())

    def test_pickle(self):
        ag_pickle = pickle.loads(pickle.dumps(self.ag))
        self.assertEqual(ag_pickle.automorphisms(), self.ag.automorphisms())


class ArchGraphSystemBugFixTest(unittest.TestCase):
    def test_duplicate_channels(self):
        def make_ag(processors, directed):
            assert processors % 2 == 0

            ag = mp.ArchGraph(directed)

            ag.add_processors(processors, 'p')

            for pe in sample(range(ag.num_processors()),
                             ag.num_processors() // 2):
                for _ in range(2):
                    ag.add_channel(pe, pe, 'L1')
                    ag.add_channel(pe, pe, 'L2')

            for pe in sample(range(ag.num_processors() - 1),
                             ag.num_processors() // 2):
                for _ in range(2):
                    ag.add_channel(pe, pe + 1, 'c')
                    ag.add_channel(pe + 1, pe, 'c')

            return ag

        ag1 = make_ag(10, True)
        self.assertEqual(ag1.num_channels(), 2 * ag1.num_processors())
        self.assertCountEqual(ag1.processor_types(), ['p', 'p%L1%L2'])
        self.assertCountEqual(ag1.channel_types(), ['c', 'L1', 'L2'])

        ag2 = make_ag(10, False)
        self.assertEqual(ag2.num_channels(), int(1.5 * ag2.num_processors()))
        self.assertCountEqual(ag2.processor_types(), ['p', 'p%L1%L2'])
        self.assertCountEqual(ag2.channel_types(), ['c', 'L1', 'L2'])

    def test_self_channels(self):
        def make_ag(directed, processors):
            assert processors % 2 == 0

            ag = mp.ArchGraph(directed)

            for _ in range(processors):
                ag.add_processor('p')

            ag.fully_connect('RAM')

            return ag

        for directed in True, False:
            # ag1
            ag1 = make_ag(directed, 10)

            for pe in range(ag1.num_processors()):
                ag1.add_channel(pe, pe, 'L1')

            self.assertEqual(ag1.num_automorphisms(),
                             factorial(ag1.num_processors()))

            # ag2
            ag2 = make_ag(directed, 20)

            for pe in range(ag2.num_processors()):
                if pe < ag2.num_processors() // 2:
                    ag2.add_channel(pe, pe, 'L1')
                else:
                    ag2.add_channel(pe, pe, 'L2')

            self.assertEqual(ag2.num_automorphisms(),
                             factorial(ag2.num_processors() // 2)**2)

            # ag 3
            ag3 = make_ag(directed, 30)

            for pe, caches in zip(range(ag3.num_processors()),
                                  cycle(permutations(['L1', 'L2', 'L3']))):
                for cache in caches:
                    ag3.add_channel(pe, pe, cache)

            self.assertEqual(ag3.num_automorphisms(),
                             factorial(ag3.num_processors()))

    def test_multi_channel_automorphisms(self):
        def make_ag(directed, processors):
            assert processors % 2 == 0

            ag = mp.ArchGraph(directed)
            ag.add_processors(processors, 'p')
            ag.fully_connect('RAM')

            return ag

        ag1 = make_ag(True, 4)
        self.assertEqual(ag1.num_automorphisms(),
                         factorial(ag1.num_processors()))

        ag2 = make_ag(False, 4)
        self.assertEqual(ag2.num_automorphisms(),
                         factorial(ag2.num_processors()))

        ag3 = make_ag(True, 4)
        ag3.add_channel(0, 1, 'c')
        ag3.add_channel(1, 2, 'c')
        ag3.add_channel(2, 3, 'c')
        ag3.add_channel(3, 0, 'c')
        self.assertEqual(ag3.num_automorphisms(), 4)

        ag4 = make_ag(True, 4)
        ag4.add_channel(0, 1, 'c')
        ag4.add_channel(1, 0, 'c')
        ag4.add_channel(1, 2, 'c')
        ag4.add_channel(2, 1, 'c')
        ag4.add_channel(2, 3, 'c')
        ag4.add_channel(3, 2, 'c')
        ag4.add_channel(3, 0, 'c')
        ag4.add_channel(0, 3, 'c')
        self.assertEqual(ag4.num_automorphisms(), 8)

        ag5 = make_ag(False, 4)
        ag5.add_channel(0, 1, 'c')
        ag5.add_channel(1, 2, 'c')
        ag5.add_channel(2, 3, 'c')
        ag5.add_channel(3, 0, 'c')
        self.assertEqual(ag5.num_automorphisms(), 8)


if __name__ == '__main__':
    unittest.main()
