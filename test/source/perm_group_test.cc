#include <vector>

#include "gmock/gmock.h"
#include "perm.h"
#include "perm_group.h"
#include "test_utility.h"

#include "test_main.cc"

using testing::ElementsAre;
using testing::UnorderedElementsAre;
using testing::UnorderedElementsAreArray;

static unsigned factorial(unsigned x)
{
  unsigned result = 1u;
  while (x > 1u)
    result *= x--;

  return result;
}

TEST(PermGroupTest, CanCalculateOrbit)
{
  cgtl::Perm perm1(5, {{1, 2, 5}});
  cgtl::Perm perm2(5, {{1, 4}, {3, 5}});
  cgtl::PermGroup pg(5, {perm1, perm2});

  cgtl::SchreierTree st(5);
  std::vector<unsigned> orbit = cgtl::PermGroup::orbit(1, {perm1, perm2}, st);

  ASSERT_THAT(orbit, UnorderedElementsAre(1, 2, 3, 4, 5))
    << "Orbit elements correct.";

  std::vector<std::vector<unsigned>> expected_transversals {
    {1, 2, 3, 4, 5}, {2, 5, 3, 4, 1}, {3, 4, 5, 1, 2}, {4, 2, 5, 1, 3},
    {5, 1, 3, 4, 2}
  };

  for (unsigned i = 1u; i <= 5u; ++i) {
    EXPECT_TRUE(perm_equal(expected_transversals[i - 1u], st.transversal(i)))
      << "Transversal for " << i << " correct.";
  }
}

TEST(PermGroupTest, SchreierSimsWorks)
{
  std::vector<unsigned> base;
  std::vector<cgtl::Perm> generators {cgtl::Perm(5, {{1, 2, 4, 3}}),
                                      cgtl::Perm(5, {{1, 2, 5, 4}})};

  std::vector<cgtl::SchreierTree> dummy;

  cgtl::PermGroup::schreier_sims(base, generators, dummy);

  EXPECT_THAT(base, ElementsAre(1, 2))
    << "Base correct.";

  EXPECT_THAT(generators, UnorderedElementsAre(
    cgtl::Perm(5, {{1, 2, 4, 3}}), cgtl::Perm(5, {{1, 2, 5, 4}}),
    cgtl::Perm(5, {{2, 5}, {3, 4}}), cgtl::Perm(5, {{2, 3, 5, 4}})))
      << "Strong generating set correct.";
}

TEST(PermGroupTest, CanObtainDegree)
{
  cgtl::PermGroup pg(10u, {cgtl::Perm(10u)});
  EXPECT_EQ(10u, pg.degree())
    << "Permutation group degree set correctly.";
}

TEST(PermGroupTest, CanObtainOrder)
{
  for (unsigned i = 1u; i <= 10u; ++i) {
    EXPECT_EQ(factorial(i), cgtl::PermGroup::symmetric(i).order())
      << "Order set correctly for symmetric group S" << i;
  }

  for (unsigned i = 1u; i <= 10u; ++i) {
    EXPECT_EQ(i, cgtl::PermGroup::cyclic(i).order())
      << "Order set correctly for cyclic group Z" << i;
  }

  for (unsigned i = 3u; i <= 10u; ++i) {
    EXPECT_EQ(factorial(i) / 2, cgtl::PermGroup::alternating(i).order())
      << "Order set correctly for alternating group A" << i;
  }
}

TEST(PermGroupTest, CanTestMembership)
{
  cgtl::PermGroup a4 = cgtl::PermGroup::alternating(4);

  std::vector<cgtl::Perm> expected_members = {
    cgtl::Perm(4),
    cgtl::Perm(4, {{2, 3, 4}}),
    cgtl::Perm(4, {{2, 4, 3}}),
    cgtl::Perm(4, {{1, 2}, {3, 4}}),
    cgtl::Perm(4, {{1, 2, 3}}),
    cgtl::Perm(4, {{1, 2, 4}}),
    cgtl::Perm(4, {{1, 3, 2}}),
    cgtl::Perm(4, {{1, 3, 4}}),
    cgtl::Perm(4, {{1, 3}, {2, 4}}),
    cgtl::Perm(4, {{1, 4, 2}}),
    cgtl::Perm(4, {{1, 4, 3}}),
    cgtl::Perm(4, {{1, 4}, {2, 3}})
  };

  std::vector<cgtl::Perm> expected_non_members = {
    cgtl::Perm(4, {{3, 4}}),
    cgtl::Perm(4, {{2, 3}}),
    cgtl::Perm(4, {{2, 4}}),
    cgtl::Perm(4, {{1, 2}}),
    cgtl::Perm(4, {{1, 2, 3, 4}}),
    cgtl::Perm(4, {{1, 2, 4, 3}}),
    cgtl::Perm(4, {{1, 3, 4, 2}}),
    cgtl::Perm(4, {{1, 3}}),
    cgtl::Perm(4, {{1, 3, 2, 4}}),
    cgtl::Perm(4, {{1, 4, 3, 2}}),
    cgtl::Perm(4, {{1, 4}}),
    cgtl::Perm(4, {{1, 4, 2, 3}})
  };

  for (cgtl::Perm const &perm : expected_members) {
    EXPECT_TRUE(a4.is_element(perm))
      << "Membership test correctly identifies group member " << perm;
  }

  for (cgtl::Perm const &perm : expected_non_members) {
    EXPECT_FALSE(a4.is_element(perm))
      << "Membership test correctly rejects non group member " << perm;
  }
}

TEST(PermGroupTest, CanGenerateRandomElement)
{
  cgtl::PermGroup a4 = cgtl::PermGroup::alternating(4);

  for (unsigned i = 0u; i < 1000u; ++i) {
    EXPECT_TRUE(a4.is_element(a4.random_element()))
      << "Randomly generated group element is actually inside group.";
  }
}

TEST(PermGroupTest, CanIterateElements)
{
  cgtl::PermGroup a4 = cgtl::PermGroup::alternating(4);

  std::vector<cgtl::Perm> expected_members = {
    cgtl::Perm(4),
    cgtl::Perm(4, {{2, 3, 4}}),
    cgtl::Perm(4, {{2, 4, 3}}),
    cgtl::Perm(4, {{1, 2}, {3, 4}}),
    cgtl::Perm(4, {{1, 2, 3}}),
    cgtl::Perm(4, {{1, 2, 4}}),
    cgtl::Perm(4, {{1, 3, 2}}),
    cgtl::Perm(4, {{1, 3, 4}}),
    cgtl::Perm(4, {{1, 3}, {2, 4}}),
    cgtl::Perm(4, {{1, 4, 2}}),
    cgtl::Perm(4, {{1, 4, 3}}),
    cgtl::Perm(4, {{1, 4}, {2, 3}})
  };

  std::vector<cgtl::Perm> actual_members1;

  for (cgtl::Perm const &p : a4)
    actual_members1.push_back(p);

  EXPECT_THAT(actual_members1, UnorderedElementsAreArray(expected_members))
    << "Iteration produces every element exactly once (ranged for).";

  std::vector<cgtl::Perm> actual_members2;

  for (cgtl::PermGroup::const_iterator it = a4.begin(); it != a4.end(); it++) {
    EXPECT_EQ(4u, it->degree())
      << "Iterator dereferencing works correctly.";

    EXPECT_TRUE(it == it && it != a4.end())
      << "Iterator comparison works correctly.";

    actual_members2.push_back(*it);
  }

  EXPECT_THAT(actual_members2, UnorderedElementsAreArray(expected_members))
    << "Iteration produces every element exactly once (explicit iterator).";
}
