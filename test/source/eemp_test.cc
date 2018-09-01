#include <utility>
#include <vector>

#include "gmock/gmock.h"

#include "eemp.h"

#include "test_main.cc"

using cgtl::EEMP;
using cgtl::PartialPerm;

using testing::ElementsAreArray;

TEST(EEMPTest, CanComputeActionComponents)
{
  std::vector<PartialPerm> generators {
    PartialPerm({4, 6, 8, 1, 5, 2, 7, 3, 9}),
    PartialPerm({5, 7, 9, 2, 4, 1, 6, 3, 8}),
    PartialPerm({0, 5, 0, 0, 6, 2}),
    PartialPerm({3, 1, 2})
  };

  EEMP::SchreierTree schreier_tree;
  EEMP::OrbitGraph orbit_graph;

  std::vector<std::vector<unsigned>> action_component =
    EEMP::action_components({1, 2, 3, 4, 5, 6, 7, 8, 9}, generators,
                            schreier_tree, orbit_graph);

  decltype(action_component) expected_action_components {
	{1, 2, 3, 4, 5, 6, 7, 8, 9}, {2, 5, 6}, {1, 2, 3}, {1, 4, 7}, {1},
    {4, 6, 8}, {5, 7, 9}, {5}, {}, {3}, {4}, {2}, {6}, {8}, {9}, {7}
  };

  decltype(schreier_tree.data) expected_schreier_tree {
	{0, 2}, {0, 3}, {1, 1}, {1, 3}, {2, 0}, {2, 1}, {2, 2}, {3, 2}, {3, 3},
    {4, 0}, {5, 2}, {6, 2}, {9, 0}, {9, 1}, {11, 1}
  };

  decltype(orbit_graph.data) expected_orbit_graph {
    {0, 1, 5, 3, 10, 2, 6, 7, 8, 13, 4, 12, 11, 9, 14, 15},
    {0, 3, 6, 1, 7, 2, 5, 10, 8, 14, 11, 15, 4, 9, 13, 12},
    {1, 1, 7, 8, 8, 11, 12, 12, 8, 8, 8, 7, 11, 8, 8, 8},
    {2, 4, 2, 9, 9, 8, 8, 8, 8, 11, 8, 4, 8, 8, 8, 8}
  };

  ASSERT_THAT(action_component, ElementsAreArray(expected_action_components))
    << "Component of action determined correctly.";

  EXPECT_THAT(schreier_tree.data, ElementsAreArray(expected_schreier_tree))
    << "Schreier tree representation correct.";

  EXPECT_THAT(orbit_graph.data, ElementsAreArray(expected_orbit_graph))
    << "Orbit graph representation correct.";

  // TODO: test schreier tree tracing
}
