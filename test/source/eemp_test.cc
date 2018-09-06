#include <sstream>
#include <utility>
#include <vector>

#include "gmock/gmock.h"

#include "eemp.h"
#include "partial_perm.h"
#include "perm.h"
#include "perm_group.h"
#include "test_utility.h"
#include "util.h"

#include "test_main.cc"

using cgtl::EEMP;
using cgtl::PartialPerm;
using cgtl::Perm;
using cgtl::PermGroup;
using cgtl::expand_partition;

using testing::ElementsAreArray;
using testing::UnorderedElementsAreArray;

class EEMPTest : public testing::Test
{
protected:
  void SetUp() {
    s = EEMP::S(
      {
        PartialPerm({4, 6, 8, 1, 5, 2, 7, 3, 9}),
        PartialPerm({5, 7, 9, 2, 4, 1, 6, 3, 8}),
        PartialPerm({0, 5, 0, 0, 6, 2}),
        PartialPerm({3, 1, 2})
      }
    );
  }

  EEMP::S s;
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

  ASSERT_THAT(s.action_component, ElementsAreArray(expected_action_component))
    << "Component of action determined correctly.";

  EXPECT_THAT(s.schreier_tree.data, ElementsAreArray(expected_schreier_tree))
    << "Schreier tree representation correct.";

  EXPECT_THAT(s.orbit_graph.data, ElementsAreArray(expected_orbit_graph))
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

  PartialPerm const x(s.generators[0] * s.generators[2] * s.generators[3]);

  std::vector<PartialPerm> inverted_generators(s.generators.size());
  for (auto i = 0u; i < s.generators.size(); ++i)
    inverted_generators[i] = ~s.generators[i];

  EEMP::SchreierTree left_schreier_tree;
  EEMP::OrbitGraph dummy;
  auto left_action_component(EEMP::action_component(
    x.dom(), inverted_generators, s.dom, left_schreier_tree, s.orbit_graph));

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

  EXPECT_THAT(expand_partition(s.scc), ElementsAreArray(expected_scc))
    << "Strongly connected components of orbit graph determined correctly.";
}

TEST_F(EEMPTest, CanTraceSchreierTree)
{
  PartialPerm const expected_pperms[] = {
    PartialPerm(9),
    s.generators[2],
    s.generators[3],
    s.generators[2] * s.generators[3],
    s.generators[2] * s.generators[1] * s.generators[2]
  };

  auto const scc(expand_partition(s.scc));

  ASSERT_EQ(sizeof(expected_pperms) / sizeof(expected_pperms[0]), scc.size())
    << "Correct number of strongly connected components.";

  for (auto i = 0u; i < scc.size(); ++i) {
    auto c_idx = scc[i][0];
    auto c(s.action_component[c_idx]);

    auto pperm(EEMP::schreier_trace(c_idx, s.schreier_tree, s.generators, s.dom));

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
  PartialPerm const pperms[] = {
    s.generators[0],
    s.generators[0] * s.generators[2],
    s.generators[0] * s.generators[3],
    s.generators[0] * s.generators[2] * s.generators[3],
    s.generators[0] * s.generators[2] * s.generators[1] * s.generators[2]
  };

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
        Perm({1})
      }
    ),
    PermGroup()
  };

  for (auto i = 0u; i < sizeof(pperms) / sizeof(pperms[0]); ++i) {
    EEMP::SchreierTree st;
    EEMP::OrbitGraph og;
    auto action_component(EEMP::action_component(
      pperms[i].im(), s.generators, s.dom, st, og));

    auto sg(EEMP::schreier_generators(
      pperms[i], s.generators, s.dom, action_component, st, og));

    EXPECT_TRUE(perm_group_equal(expected_groups[i], sg))
      << "Obtained correct schreier generator generating set.";
  }
}

TEST_F(EEMPTest, CanObtainRClassElements)
{
  auto x(s.generators[0] * s.generators[2] * s.generators[3]);
  auto r(EEMP::r_class_elements(x, s.generators, s.dom));

  PartialPerm const expected_r_class[] = {
    PartialPerm({2}, {1}),
    PartialPerm({2}, {2}),
    PartialPerm({2}, {3}),
    PartialPerm({2}, {4}),
    PartialPerm({2}, {5}),
    PartialPerm({2}, {6}),
    PartialPerm({2}, {7}),
    PartialPerm({2}, {8}),
    PartialPerm({2}, {9})
  };

  EXPECT_THAT(r, UnorderedElementsAreArray(expected_r_class))
    << "Computed R-class elements are correct.";
}

TEST_F(EEMPTest, CanObtainRClassesInDClass)
{
  auto x(s.generators[0] * s.generators[2] * s.generators[3]);
  auto r(EEMP::r_classes_in_d_class(x, s.generators, s.dom));

  PartialPerm const expected_r_class_repr[] = {
    PartialPerm({2}, {1}),
    PartialPerm({3}, {1}),
    PartialPerm({1}, {1}),
    PartialPerm({6}, {1}),
    PartialPerm({7}, {1}),
    PartialPerm({8}, {1}),
    PartialPerm({4}, {1}),
    PartialPerm({5}, {1}),
    PartialPerm({9}, {1})
  };

  EXPECT_THAT(r, UnorderedElementsAreArray(expected_r_class_repr))
    << "Computed R-class representatives are correct.";
}

TEST_F(EEMPTest, CanTestRShipMembership)
{
  PartialPerm y({1, 2, 3}, {5, 7, 9});

  std::pair<PartialPerm, bool> const expected_elem_tests[] = {
    {s.generators[2], false},
    {s.generators[3], true}
  };

  for (auto const &elem_test : expected_elem_tests) {
    PartialPerm const &x(std::get<0>(elem_test));
    EEMP::RClass r_class(x, s);

    EXPECT_EQ(std::get<1>(elem_test), r_class.is_member(y))
      << "Can determine R class membership.";
  }
}
