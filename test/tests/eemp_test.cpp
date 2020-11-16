#if 0
#include <sstream>
#include <utility>
#include <vector>

#include "gmock/gmock.h"

#include "eemp.hpp"
#include "partial_perm.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "test_utility.hpp"

#include "test_main.cpp"

using namespace mpsym;
using namespace mpsym::internal;
using namespace mpsym::internal::eemp;

using testing::ElementsAreArray;
using testing::UnorderedElementsAreArray;

class EEMPTest : public testing::Test
{
protected:
  void SetUp() {
    component = action_component(
      dom, gens, dom.back(), schreier_tree, orbit_graph);

    auto tmp(strongly_connected_components(orbit_graph));
    scc = tmp.second;
    scc_expanded = expand_partition(scc);
  }

  std::vector<unsigned> const dom {1, 2, 3, 4, 5, 6, 7, 8, 9};

  std::vector<PartialPerm> const gens {
    PartialPerm({4, 6, 8, 1, 5, 2, 7, 3, 9}),
    PartialPerm({5, 7, 9, 2, 4, 1, 6, 3, 8}),
    PartialPerm({0, 5, 0, 0, 6, 2}),
    PartialPerm({3, 1, 2})
  };

  std::vector<PartialPerm> const inv_gens {
    ~PartialPerm({4, 6, 8, 1, 5, 2, 7, 3, 9}),
    ~PartialPerm({5, 7, 9, 2, 4, 1, 6, 3, 8}),
    ~PartialPerm({0, 5, 0, 0, 6, 2}),
    ~PartialPerm({3, 1, 2})
  };

  std::vector<std::vector<unsigned>> component;
  eemp::SchreierTree schreier_tree;
  OrbitGraph orbit_graph;
  std::vector<unsigned> scc;
  std::vector<std::vector<unsigned>> scc_expanded;
};

TEST_F(EEMPTest, CanComputeActionComponent)
{
  std::vector<unsigned> const expected_action_component[] = {
	{1, 2, 3, 4, 5, 6, 7, 8, 9}, {2, 5, 6}, {1, 2, 3}, {1, 4, 7}, {1},
    {4, 6, 8}, {5, 7, 9}, {5}, {}, {3}, {4}, {2}, {6}, {8}, {9}, {7}
  };

  std::pair<unsigned, unsigned> const expected_schreier_tree[] = {
	{0, 2}, {0, 3}, {1, 1}, {1, 3}, {2, 0}, {2, 1}, {2, 2}, {3, 2}, {3, 3},
    {4, 0}, {5, 2}, {6, 2}, {9, 0}, {9, 1}, {11, 1}
  };

  std::vector<unsigned> const expected_orbit_graph[] = {
    {0, 1, 5, 3, 10, 2, 6, 7, 8, 13, 4, 12, 11, 9, 14, 15},
    {0, 3, 6, 1, 7, 2, 5, 10, 8, 14, 11, 15, 4, 9, 13, 12},
    {1, 1, 7, 8, 8, 11, 12, 12, 8, 8, 8, 7, 11, 8, 8, 8},
    {2, 4, 2, 9, 9, 8, 8, 8, 8, 11, 8, 4, 8, 8, 8, 8}
  };

  ASSERT_THAT(component, ElementsAreArray(expected_action_component))
    << "Component of action determined correctly.";

  EXPECT_THAT(schreier_tree.data, ElementsAreArray(expected_schreier_tree))
    << "Schreier tree representation correct.";

  EXPECT_THAT(orbit_graph.data, ElementsAreArray(expected_orbit_graph))
    << "Orbit graph representation correct.";
}

TEST_F(EEMPTest, CanComputeLeftSchreierTree)
{
  std::vector<unsigned> const expected_left_action_component[] = {
	{2}, {6}, {4}, {3}, {7}, {5}, {}, {1}, {8}, {9}
  };

  std::pair<unsigned, unsigned> const expected_left_schreier_tree[] = {
	{0, 0}, {0, 1}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {3, 0}, {8, 1}
  };

  PartialPerm const x(gens[0] * gens[2] * gens[3]);

  eemp::SchreierTree left_schreier_tree;
  OrbitGraph dummy;
  auto left_action_component(action_component(
    x.dom(), inv_gens, dom.back(), left_schreier_tree, orbit_graph));

  ASSERT_THAT(left_action_component,
              ElementsAreArray(expected_left_action_component))
    << "Component of action determined correctly.";

  EXPECT_THAT(left_schreier_tree.data,
              ElementsAreArray(expected_left_schreier_tree))
    << "Schreier tree representation correct.";
}

TEST_F(EEMPTest, CanIdentifyStronglyConnectedOrbitGraphComponents)
{
  std::vector<unsigned> const expected_scc[] = {
    {0}, {1, 3}, {2, 5, 6}, {4, 7, 9, 10, 11, 12, 13, 14, 15}, {8}
  };

  EXPECT_THAT(scc_expanded, ElementsAreArray(expected_scc))
    << "Strongly connected components of orbit graph determined correctly.";
}

TEST_F(EEMPTest, CanComputeSccSpanningTrees)
{
  scc_spanning_tree(0, orbit_graph, scc);
  scc_spanning_tree(1, orbit_graph, scc);
  scc_spanning_tree(2, orbit_graph, scc);
  scc_spanning_tree(4, orbit_graph, scc);
}

TEST_F(EEMPTest, CanTraceSchreierTree)
{
  PartialPerm const expected_pperms[] = {
    PartialPerm(9),
    gens[2],
    gens[3],
    gens[2] * gens[3],
    gens[2] * gens[1] * gens[2]
  };

  for (auto i = 0u; i < scc_expanded.size(); ++i) {
    auto c_idx = scc_expanded[i][0];
    auto c(component[c_idx]);

    auto pperm(schreier_trace(c_idx, schreier_tree, gens, dom.back()));

    std::stringstream ss;

    ss << '[';
    for (auto j = 0u; j < c.size(); ++j) {
      ss << c[j];
      if (j != c.size() - 1u)
        ss << ", ";
    }
    ss << ']';

    EXPECT_EQ(expected_pperms[i], pperm)
      << "Partial permutation determined for action component "
      << c_idx + 1u << " (" << ss.str() << ") traced correctly.";
  }
}

TEST_F(EEMPTest, CanComputeStabilizerSchreierGenerators)
{
  unsigned const scc_repr[] = {0, 1, 2, 4, 8};

  PermGroup const expected_groups[] = {
    PermGroup(9,
      {
        Perm(9, {{1, 4}, {2, 6}, {3, 8}}),
        Perm(9, {{1, 5, 4, 2, 7, 6}, {3, 9, 8}})
      }
    ),
    PermGroup(6,
      {
        Perm(6, {{2, 6}}),
        Perm(6, {{2, 6, 5}})
      }
    ),
    PermGroup(3,
      {
        Perm(3, {{1, 3, 2}}),
        Perm(3, {{1, 2}})
      }
    ),
    PermGroup(1,
      {
        Perm{1}
      }
    ),
    PermGroup()
  };

  for (auto i = 0u; i < sizeof(scc_repr) / sizeof(scc_repr[0]); ++i) {
    eemp::SchreierTree st(scc_spanning_tree(scc_repr[i], orbit_graph, scc));

    auto sg(schreier_generators(
      scc_repr[i], gens, dom.back(), component, st, orbit_graph, scc));

    EXPECT_TRUE(perm_group_equal(expected_groups[i], sg))
      << "Obtained correct schreier generator generating set.";
  }
}

TEST_F(EEMPTest, CanObtainRClassRepresentatives)
{
  PartialPerm const expected_r_class_repr[] = {
    gens[0],
    gens[0] * gens[2],
    gens[0] * gens[3],
    gens[0] * gens[2] * gens[1],
    gens[0] * gens[2] * gens[3],
    gens[0] * gens[3] * gens[0],
    gens[0] * gens[3] * gens[1],
    gens[0] * gens[3] * gens[2],
    gens[0] * gens[2] * gens[1] * gens[2],
    gens[0] * gens[2] * gens[1] * gens[3],
    gens[0] * gens[2] * gens[3] * gens[0],
    gens[0] * gens[3] * gens[0] * gens[2],
    gens[0] * gens[3] * gens[1] * gens[2],
    gens[0] * gens[2] * gens[1] * gens[3] * gens[0],
    gens[0] * gens[2] * gens[1] * gens[3] * gens[1],
    gens[0] * gens[3] * gens[0] * gens[2] * gens[1]
  };

  auto r_class_repr(r_class_representatives(schreier_tree, gens));

  EXPECT_THAT(r_class_repr, UnorderedElementsAreArray(expected_r_class_repr))
    << "R class representatives determined correctly.";
}
#endif
