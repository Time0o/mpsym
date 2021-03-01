import pickle
import unittest
from copy import deepcopy
from itertools import cycle, permutations
from math import factorial
from random import sample
from textwrap import dedent

import mpsym as mp


def cyclically_connect(self, ct, start=0, end=None, directed=True):
    if end is None:
        end = self.num_processors() - 1

    def connect(source, target):
        self.add_channel(source, target, ct)

        if self.directed() and not directed:
            self.add_channel(target, source, ct)

    for pe in range(start, end):
        connect(pe, pe + 1)

    connect(end, start)

mp.ArchGraph.cyclically_connect = cyclically_connect


class PermTest(unittest.TestCase):
    def setUp(self):
        self.dom = [0,1,2,3]

        self.p = mp.Perm([1,0,3,2])
        self.q = mp.Perm([1,2,0,3])

    def test_equality(self):
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
        self.assertTrue(self.ag.num_automorphisms() == len(self.ag.automorphisms()) == 8192)

    def test_representative(self):
        for orbit in [self.ag_orbit1, self.ag_orbit2]:
            for mapping in orbit:
                for method in 'iterate', 'orbit':
                    self.assertEqual(self.ag.representative(mapping, method=method), orbit[0])

    def test_orbit(self):
        for orbit in [self.ag_orbit1, self.ag_orbit2]:
            self.assertCountEqual(list(self.ag.orbit(orbit[0])), orbit)

        def orbit_len(orb):
            return sum(1 for _ in orb)

        for n in range(3, 7):
          Sn = mp.PermGroup.symmetric(n)
          ag = mp.ArchGraphAutomorphisms(Sn)

          self.assertEqual(orbit_len(ag.orbit(range(n))), factorial(n))

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
        self.assertEqual(ag1.processor_types(), ['p'])
        self.assertEqual(ag1.channel_types(), ['L1', 'L2', 'c'])

        ag2 = make_ag(10, False)
        self.assertEqual(ag2.num_channels(), int(1.5 * ag2.num_processors()))
        self.assertCountEqual(ag2.processor_types(), ['p'])
        self.assertCountEqual(ag2.channel_types(), ['L1', 'L2', 'c'])

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

    def test_architectures(self):
        def make_ag(processors, directed=True, mult='single'):
            assert processors % 2 == 0

            def make_ag_():
                ag = mp.ArchGraph(directed)

                if mult == 'single':
                    ag.add_processors(processors, 'p')
                elif mult == 'double_identical':
                    ag.add_processors(processors, 'p')
                    ag.add_processors(processors, 'p')
                elif mult == 'double_distinct':
                    ag.add_processors(processors, 'p1')
                    ag.add_processors(processors, 'p2')
                else:
                    raise ValueError

                return ag

            return make_ag_

        def ram(ag):
            ag.fully_connect('RAM')

        def cache(where, cache_type):
            return lambda ag: ag.fully_connect(where, cache_type)

        def ring(directed=False, partial=False):
            def ring_(ag):
                end = (ag.num_processors() // 2 - 1) if partial else None

                ag.cyclically_connect('c', start=0, end=end, directed=directed)

            return ring_

        self._test_ag_generator(make_ag(4),
                                [ram],
                                factorial(4))

        self._test_ag_generator(make_ag(4, directed=False),
                                [ram],
                                factorial(4))

        self._test_ag_generator(make_ag(4),
                                [ram, ring(directed=True)],
                                4)

        self._test_ag_generator(make_ag(4),
                                [ram, ring()],
                                8)

        self._test_ag_generator(make_ag(4, directed=False),
                                [ram, ring()],
                                8)

        self._test_ag_generator(make_ag(8),
                                [ram, ring(directed=True, partial=True)],
                                4 * factorial(4))

        self._test_ag_generator(make_ag(8),
                                [ram, ring(partial=True)],
                                8 * factorial(4))

        self._test_ag_generator(make_ag(8, directed=False),
                                [ram, ring(partial=True)],
                                8 * factorial(4))

        self._test_ag_generator(make_ag(4, mult='double_identical'),
                                [ram],
                                factorial(8))

        self._test_ag_generator(make_ag(4, mult='double_identical'),
                                [ram,
                                 cache(range(4), 'L1_1'),
                                 cache(range(4, 8), 'L1_2')],
                                factorial(4)**2)

        self._test_ag_generator(make_ag(4, mult='double_distinct'),
                                [ram],
                                factorial(4)**2)

        self._test_ag_generator(make_ag(4, mult='double_distinct'),
                                [ram,
                                 cache('p1', 'L1_1'),
                                 cache('p2', 'L1_2')],
                                factorial(4)**2)

        self._test_ag_generator(make_ag(4, mult='single'),
                                [ram,
                                 cache('p', 'L1'),
                                 ring()],
                                8)

        self._test_ag_generator(make_ag(4, mult='double_distinct'),
                                [ram,
                                 cache('p1', 'L1_1'),
                                 ring(partial=True)],
                                8 * factorial(4))

    def test_super_graph_architectures(self):
        # trivial proto or super graph
        ag1 = mp.ArchGraph(directed=False)
        ag1.add_processor('p1')
        ag1.add_processor('p2')

        ag2 = mp.ArchGraph(directed=False)
        ag2.add_processors(2, 'p1')
        ag2.fully_connect('c')

        ag3 = mp.ArchGraph(directed=False)
        ag3.add_processors(2, 'p1')
        ag3.add_processor('p2')
        ag3.fully_connect('c')

        self._test_super_graph(ag1, ag1, ['()'])
        self._test_super_graph(ag1, ag2, ['(0,1)(2,3)'])
        self._test_super_graph(ag1, ag3, ['(0,1)(3,4)'])
        self._test_super_graph(ag2, ag1, ['(0,2)(1,3)'])
        self._test_super_graph(ag2, ag2, ['(0,2)(1,3)', '(0,1)'])
        self._test_super_graph(ag2, ag3, ['(0,3)(1,4)(2,5)','(3,4)'])
        self._test_super_graph(ag3, ag1, ['(0,2)(1,3)'])
        self._test_super_graph(ag3, ag2, ['(0,2)(1,3)', '(0,1)', '(4,5)'])
        self._test_super_graph(ag3, ag3, ['(0,3)(1,4)(2,5)', '(0,1)', '(6,7)'])

    def test_symmetric_graphs(self):
        ag = mp.ArchGraph()
        ag.add_processors(10, 'p')
        ag.fully_connect('c')

        self.assertEqual(ag.num_automorphisms(), factorial(10))

        representatives = [
            ((0, 1, 2, 3), (0, 1, 2, 3)),
            ((1, 2, 3, 0), (0, 1, 2, 3)),
            ((2, 3, 0, 1), (0, 1, 2, 3)),
            ((3, 0, 1, 2), (0, 1, 2, 3)),
            ((1, 1, 2, 2), (0, 0, 1, 1)),
            ((2, 1, 2, 1), (0, 1, 0, 1)),
            ((1, 2, 3, 3), (0, 1, 2, 2)),
            ((3, 3, 3, 3), (0, 0, 0, 0)),
            ((1, 1, 10, 10), (0, 0, 10, 10)),
            ((10, 1, 10, 1), (10, 0, 10, 0)),
            ((1, 2, 10, 10), (0, 1, 10, 10)),
            ((10, 10, 10, 10), (10, 10, 10, 10))
        ]

        for mapping, representative in representatives:
            self.assertEqual(ag.representative(mapping), representative)

        ag = mp.ArchGraph()
        ag.add_processors(3, 'p1')
        ag.add_processor('p2')
        ag.add_processors(3, 'p1')
        ag.add_processor('p3')
        ag.add_processors(4, 'p1')
        ag.fully_connect('p1', 'c')

        self.assertEqual(ag.num_automorphisms(), factorial(10))

        representatives = [
            ((0, 1, 2, 3), (0, 1, 2, 3)),
            ((1, 2, 3, 0), (0, 1, 3, 2)),
            ((2, 3, 0, 1), (0, 3, 1, 2)),
            ((3, 0, 1, 2), (3, 0, 1, 2)),
            ((1, 1, 2, 2), (0, 0, 1, 1)),
            ((2, 1, 2, 1), (0, 1, 0, 1)),
            ((1, 2, 3, 3), (0, 1, 3, 3)),
            ((1, 1, 3, 7), (0, 0, 3, 7)),
            ((3, 1, 7, 1), (3, 0, 7, 0)),
            ((1, 2, 3, 7), (0, 1, 3, 7)),
            ((3, 7, 3, 7), (3, 7, 3, 7))
        ]

        for mapping, representative in representatives:
            self.assertEqual(ag.representative(mapping), representative)

    def test_kalray_automorphisms(self):
        def kalray_like_lua(naive,
                            num_processors,
                            num_clusters,
                            cross_connected=False,
                            directed=False):

            shared = dedent(
                f"""
                local mpsym = require 'mpsym'

                local NUM_PROCESSORS = {num_processors}
                local NUM_CLUSTERS = {num_clusters}

                local CROSS_CONNECTED = {'true' if cross_connected else 'false'}

                local DIRECTED = {'true' if directed else 'false'}
                """
            )

            if naive:
                return shared + dedent(
                    """
                    local processors = {}
                    local channels = {}

                    local socs = {}

                    for i = 1,NUM_CLUSTERS do
                      local cluster_processors = mpsym.combine_processors(
                        mpsym.identical_processors(NUM_PROCESSORS, 'P'), {NUM_PROCESSORS + 1, 'PS'})

                      socs[i] = mpsym.append_processors(processors, cluster_processors)
                    end

                    for i = 1,NUM_CLUSTERS do
                      mpsym.append_channels(channels, mpsym.fully_connected_channels(socs[i], 'shared memory'))
                    end

                    if CROSS_CONNECTED then
                        for i = 1,NUM_CLUSTERS do
                          for j = i+1,NUM_CLUSTERS do
                            mpsym.append_channels(channels, mpsym.fully_cross_connected_channels(socs[i], socs[j], 'C'))
                          end
                        end
                    else
                        mpsym.append_channels(channels, mpsym.fully_connected_channels(processors, 'C'))
                    end

                    return mpsym.ArchGraph:create{
                      directed = DIRECTED,
                      processors = processors,
                      channels = channels
                    }
                    """
                )
            else:
                return shared + dedent(
                    """
                    local super_graph_clusters = mpsym.identical_clusters(NUM_CLUSTERS, 'SoC')
                    local super_graph_channels = mpsym.fully_connected_channels(super_graph_clusters, 'C')

                    local proto_processors = mpsym.combine_processors(
                        mpsym.identical_processors(NUM_PROCESSORS, 'P'), {NUM_PROCESSORS + 1, 'PS'})
                    local proto_channels = mpsym.fully_connected_channels(proto_processors, 'shared memory')

                    return mpsym.ArchUniformSuperGraph:create{
                      super_graph = mpsym.ArchGraph:create{
                        directed = false,
                        clusters = super_graph_clusters,
                        channels = super_graph_channels
                      },
                      proto = mpsym.ArchGraph:create{
                        directed = DIRECTED,
                        processors = proto_processors,
                        channels = proto_channels
                      }
                    }
                    """
                )

        def kalray_like(**kwargs):
            return mp.ArchGraphSystem.from_lua(
                kalray_like_lua(num_processors=8, num_clusters=5, **kwargs))

        kalray_like_variants = [
            {'naive': False, 'directed': False},
            {'naive': False, 'directed': True},
            {'naive': True, 'cross_connected': True, 'directed': False},
            {'naive': True, 'cross_connected': False, 'directed': False},
            {'naive': True, 'cross_connected': True, 'directed': True},
            {'naive': True, 'cross_connected': False, 'directed': True}
        ]

        kalray_num_automorphisms = set()

        for kwargs in kalray_like_variants:
            kalray = kalray_like(**kwargs)
            kalray_num_automorphisms.add(kalray.num_automorphisms())

        self.assertEqual(len(kalray_num_automorphisms), 1)

    def _test_ag_generator(self, make_ag, connections, num_automorphisms):
        for fs in permutations(connections):
            ag = make_ag()

            p = 1

            for ag_ in ag, deepcopy(ag):
                for f in fs:
                    f(ag_)

                self.assertEqual(ag_.num_automorphisms(),
                                 num_automorphisms,
                                 "(permutation {})".format(p))

                p += 1

    def _test_super_graph(self, super_graph, proto, automorphism_generators):
        assert not super_graph.directed()
        assert not proto.directed()

        # graph properties
        ag = mp.ArchUniformSuperGraph(super_graph, proto)

        self.assertEqual(ag.num_processors(),
                         proto.num_processors() * super_graph.num_processors())

        self.assertEqual(ag.num_channels(),
                         proto.num_channels() * super_graph.num_processors() + \
                         super_graph.num_channels() * (proto.num_processors()**2))

        # automorphism properties
        ag_automs = ag.automorphisms()

        expected_automs = mp.PermGroup([mp.Perm(ag.num_processors(), gen)
                                        for gen in automorphism_generators])

        self.assertEqual(ag_automs, expected_automs)
        self.assertEqual(ag_automs.degree(), ag.num_processors())
        self.assertEqual(len(ag_automs), ag.num_automorphisms())

        # representatives
        ag_automs_graph = mp.ArchGraphAutomorphisms(ag_automs)

        id_mapping = tuple(i for i in range(min(ag.num_processors(), 5)))

        for mapping in permutations(id_mapping):
            self.assertEqual(ag.representative(mapping),
                             ag_automs_graph.representative(mapping))
