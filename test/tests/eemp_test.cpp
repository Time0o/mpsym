#include <sstream>
#include <utility>
#include <vector>

#include "gmock/gmock.h"

#include "eemp.hpp"
#include "partial_perm.hpp"
#include "partial_perm_set.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "test_utility.hpp"

#include "test_main.cpp"

using namespace mpsym;
using namespace mpsym::internal;
using namespace mpsym::internal::eemp;

using testing::ElementsAreArray;
using testing::UnorderedElementsAreArray;

namespace
{

class EEMPTest : public testing::Test
{
protected:
  void SetUp() {
    component = action_component(dom, gens, schreier_tree, orbit_graph);

    sccs = strongly_connected_components(orbit_graph);
    sccs_expanded = sccs.data_expanded();
  }

  std::vector<unsigned> const dom {0, 1, 2, 3, 4, 5, 6, 7, 8};

  PartialPermSet const gens {
    PartialPerm({3, 5, 7, 0, 4, 1, 6, 2, 8}),
    PartialPerm({4, 6, 8, 1, 3, 0, 5, 2, 7}),
    PartialPerm({-1, 4, -1, -1, 5, 1}),
    PartialPerm({2, 0, 1})
  };

  PartialPermSet const inv_gens {
    ~PartialPerm({3, 5, 7, 0, 4, 1, 6, 2, 8}),
    ~PartialPerm({4, 6, 8, 1, 3, 0, 5, 2, 7}),
    ~PartialPerm({-1, 4, -1, -1, 5, 1}),
    ~PartialPerm({2, 0, 1})
  };

  eemp::Component component;
  eemp::SchreierTree schreier_tree;
  eemp::OrbitGraph orbit_graph;

  eemp::Sccs sccs;
  std::vector<std::vector<unsigned>> sccs_expanded;
};

TEST_F(EEMPTest, CanComputeActionComponent)
{
  std::vector<unsigned> const expected_action_component[] = {
	{0, 1, 2, 3, 4, 5, 6, 7, 8}, {1, 4, 5}, {0, 1, 2}, {0, 3, 6}, {0},
    {3, 5, 7}, {4, 6, 8}, {4}, {}, {2}, {3}, {1}, {5}, {7}, {8}, {6}
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
	{1}, {5}, {3}, {2}, {6}, {4}, {}, {0}, {7}, {8}
  };

  std::pair<unsigned, unsigned> const expected_left_schreier_tree[] = {
	{0, 0}, {0, 1}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {3, 0}, {8, 1}
  };

  PartialPerm const x(gens[0] * gens[2] * gens[3]);

  eemp::SchreierTree left_schreier_tree;
  OrbitGraph dummy;
  auto left_action_component(action_component(
    x.dom(), inv_gens, left_schreier_tree, orbit_graph));

  ASSERT_THAT(left_action_component,
              ElementsAreArray(expected_left_action_component))
    << "Component of action determined correctly.";

  EXPECT_THAT(left_schreier_tree.data,
              ElementsAreArray(expected_left_schreier_tree))
    << "Schreier tree representation correct.";
}

// TODO: remove?
TEST_F(EEMPTest, CanIdentifyStronglyConnectedOrbitGraphComponents)
{
  std::stringstream ss;
  ss << sccs;

  EXPECT_EQ("{{0}, {1, 3}, {2, 5, 6}, {4, 7, 9, 10, 11, 12, 13, 14, 15}, {8}}", ss.str())
    << "Strongly connected components of orbit graph determined correctly.";
}

// TODO: remove?
TEST_F(EEMPTest, CanComputeSccSpanningTrees)
{
  spanning_tree(orbit_graph, sccs, 0);
  spanning_tree(orbit_graph, sccs, 1);
  spanning_tree(orbit_graph, sccs, 2);
  spanning_tree(orbit_graph, sccs, 4);
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

  for (auto i = 0u; i < sccs_expanded.size(); ++i) {
    auto c_idx = sccs_expanded[i][0];
    auto c(component[c_idx]);

    auto pperm(schreier_trace(gens, schreier_tree, c_idx, dom.back()));

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

//TEST_F(EEMPTest, CanComputeStabilizerSchreierGenerators)
//{
//  unsigned const scc_repr[] = {0, 1, 2, 4, 8};
//
//  PermGroup const expected_groups[] = {
//    PermGroup(9,
//      {
//        Perm(9, {{1, 4}, {2, 6}, {3, 8}}),
//        Perm(9, {{1, 5, 4, 2, 7, 6}, {3, 9, 8}})
//      }
//    ),
//    PermGroup(6,
//      {
//        Perm(6, {{2, 6}}),
//        Perm(6, {{2, 6, 5}})
//      }
//    ),
//    PermGroup(3,
//      {
//        Perm(3, {{1, 3, 2}}),
//        Perm(3, {{1, 2}})
//      }
//    ),
//    PermGroup(1,
//      {
//        Perm{1}
//      }
//    ),
//    PermGroup()
//  };
//
//  for (auto i = 0u; i < sizeof(scc_repr) / sizeof(scc_repr[0]); ++i) {
//    eemp::SchreierTree st(scc_spanning_tree(scc_repr[i], orbit_graph, scc));
//
//    auto sg(schreier_generators(
//      scc_repr[i], gens, dom.back(), component, st, orbit_graph, scc));
//
//    EXPECT_TRUE(perm_group_equal(expected_groups[i], sg))
//      << "Obtained correct schreier generator generating set.";
//  }
//}
//
//TEST_F(EEMPTest, CanObtainRClassRepresentatives)
//{
//  PartialPerm const expected_r_class_repr[] = {
//    gens[0],
//    gens[0] * gens[2],
//    gens[0] * gens[3],
//    gens[0] * gens[2] * gens[1],
//    gens[0] * gens[2] * gens[3],
//    gens[0] * gens[3] * gens[0],
//    gens[0] * gens[3] * gens[1],
//    gens[0] * gens[3] * gens[2],
//    gens[0] * gens[2] * gens[1] * gens[2],
//    gens[0] * gens[2] * gens[1] * gens[3],
//    gens[0] * gens[2] * gens[3] * gens[0],
//    gens[0] * gens[3] * gens[0] * gens[2],
//    gens[0] * gens[3] * gens[1] * gens[2],
//    gens[0] * gens[2] * gens[1] * gens[3] * gens[0],
//    gens[0] * gens[2] * gens[1] * gens[3] * gens[1],
//    gens[0] * gens[3] * gens[0] * gens[2] * gens[1]
//  };
//
//  auto r_class_repr(r_class_representatives(schreier_tree, gens));
//
//  EXPECT_THAT(r_class_repr, UnorderedElementsAreArray(expected_r_class_repr))
//    << "R class representatives determined correctly.";
//}

} // anonymous namespace
