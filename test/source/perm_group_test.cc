#include <sstream>
#include <utility>
#include <vector>

#include "gmock/gmock.h"
#include "perm.h"
#include "perm_group.h"
#include "schreier_sims.h"
#include "test_utility.h"

#include "test_main.cc"

using cgtl::Perm;
using cgtl::PermGroup;
using cgtl::SchreierSims;

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

TEST(PermGroupTest, CanObtainDegree)
{
  PermGroup pg(10u, {Perm(10u)});
  EXPECT_EQ(10u, pg.degree())
    << "Permutation group degree set correctly.";
}

TEST(PermGroupTest, CanObtainOrder)
{
  PermGroup id(10u, {});
  EXPECT_EQ(1u, id.order())
    << "Order set correctly for trivial permutation group.";

  for (unsigned i = 1u; i <= 10u; ++i) {
    EXPECT_EQ(factorial(i), PermGroup::symmetric(i).order())
      << "Order set correctly for symmetric group S" << i;
  }

  for (unsigned i = 1u; i <= 10u; ++i) {
    EXPECT_EQ(i, PermGroup::cyclic(i).order())
      << "Order set correctly for cyclic group Z" << i;
  }

  for (unsigned i = 3u; i <= 10u; ++i) {
    EXPECT_EQ(factorial(i) / 2, PermGroup::alternating(i).order())
      << "Order set correctly for alternating group A" << i;
  }
}

TEST(PermGroupTest, CanTestMembership)
{
  PermGroup a4 = PermGroup::alternating(4);

  std::vector<Perm> expected_members = {
    Perm(4),
    Perm(4, {{2, 3, 4}}),
    Perm(4, {{2, 4, 3}}),
    Perm(4, {{1, 2}, {3, 4}}),
    Perm(4, {{1, 2, 3}}),
    Perm(4, {{1, 2, 4}}),
    Perm(4, {{1, 3, 2}}),
    Perm(4, {{1, 3, 4}}),
    Perm(4, {{1, 3}, {2, 4}}),
    Perm(4, {{1, 4, 2}}),
    Perm(4, {{1, 4, 3}}),
    Perm(4, {{1, 4}, {2, 3}})
  };

  std::vector<Perm> expected_non_members = {
    Perm(4, {{3, 4}}),
    Perm(4, {{2, 3}}),
    Perm(4, {{2, 4}}),
    Perm(4, {{1, 2}}),
    Perm(4, {{1, 2, 3, 4}}),
    Perm(4, {{1, 2, 4, 3}}),
    Perm(4, {{1, 3, 4, 2}}),
    Perm(4, {{1, 3}}),
    Perm(4, {{1, 3, 2, 4}}),
    Perm(4, {{1, 4, 3, 2}}),
    Perm(4, {{1, 4}}),
    Perm(4, {{1, 4, 2, 3}})
  };

  for (Perm const &perm : expected_members) {
    EXPECT_TRUE(a4.is_element(perm))
      << "Membership test correctly identifies group member " << perm;
  }

  for (Perm const &perm : expected_non_members) {
    EXPECT_FALSE(a4.is_element(perm))
      << "Membership test correctly rejects non group member " << perm;
  }
}

TEST(PermGroupTest, CanGenerateRandomElement)
{
  PermGroup a4 = PermGroup::alternating(4);

  for (unsigned i = 0u; i < 1000u; ++i) {
    EXPECT_TRUE(a4.is_element(a4.random_element()))
      << "Randomly generated group element is actually inside group.";
  }
}

TEST(PermGroupTest, CanIterateTrivialGroup)
{
  PermGroup id = PermGroup(4, {});

  std::vector<Perm> actual_members1;
  for (Perm const &p : id)
    actual_members1.push_back(p);

  ASSERT_EQ(1u, actual_members1.size())
    << "Iterating trivial permutation group yields one element (ranged for).";

  EXPECT_TRUE(perm_equal({1, 2, 3, 4}, actual_members1[0]))
    << "Iterating trivial permutation group yields identity permutatio (ranged for)..";

  std::vector<Perm> actual_members2;
  for (PermGroup::const_iterator it = id.begin(); it != id.end(); it++)
    actual_members2.push_back(*it);

  ASSERT_EQ(1u, actual_members2.size())
    << "Iterating trivial permutation group yields one element (explicit iterator).";

  EXPECT_TRUE(perm_equal({1, 2, 3, 4}, actual_members2[0]))
    << "Iterating trivial permutation group yields identity permutation (explicit iterator).";
}

TEST(PermGroupTest, CanIterateElements)
{
  PermGroup a4 = PermGroup::alternating(4);

  std::vector<Perm> expected_members = {
    Perm(4),
    Perm(4, {{2, 3, 4}}),
    Perm(4, {{2, 4, 3}}),
    Perm(4, {{1, 2}, {3, 4}}),
    Perm(4, {{1, 2, 3}}),
    Perm(4, {{1, 2, 4}}),
    Perm(4, {{1, 3, 2}}),
    Perm(4, {{1, 3, 4}}),
    Perm(4, {{1, 3}, {2, 4}}),
    Perm(4, {{1, 4, 2}}),
    Perm(4, {{1, 4, 3}}),
    Perm(4, {{1, 4}, {2, 3}})
  };

  std::vector<Perm> actual_members1;

  for (Perm const &p : a4)
    actual_members1.push_back(p);

  EXPECT_THAT(actual_members1, UnorderedElementsAreArray(expected_members))
    << "Iteration produces every element exactly once (ranged for).";

  std::vector<Perm> actual_members2;

  for (PermGroup::const_iterator it = a4.begin(); it != a4.end(); it++) {
    EXPECT_EQ(4u, it->degree())
      << "Iterator dereferencing works correctly.";

    EXPECT_TRUE(it == it && it != a4.end())
      << "Iterator comparison works correctly.";

    actual_members2.push_back(*it);
  }

  EXPECT_THAT(actual_members2, UnorderedElementsAreArray(expected_members))
    << "Iteration produces every element exactly once (explicit iterator).";
}

class SchreierSimsVariantTest :
  public testing::TestWithParam<SchreierSims::Variant> {};

// TODO: test more groups
TEST_P(SchreierSimsVariantTest, CanGenerateCorrectGroupElements)
{
  SchreierSims::Variant schreier_var = GetParam();

  typedef std::vector<std::vector<std::vector<unsigned>>> elemset;

  std::vector<std::pair<unsigned, elemset>> groups {
    {4, {{{2, 4}}, {{1, 2}, {3, 4}}}}
  };

  std::vector<elemset> expected_elements {
    {
      {{1, 2, 3, 4}},
      {{1, 2}, {3, 4}},
      {{1, 3}, {2, 4}},
      {{1, 3}},
      {{1, 4, 3, 2}},
      {{1, 4}, {2, 3}},
      {{2, 4}}
    }
  };

  for (auto i = 0u; i < groups.size(); ++i) {
    unsigned degree = std::get<0>(groups[i]);

    std::vector<Perm> generators;
    for (auto const &perm : std::get<1>(groups[i]))
      generators.push_back(Perm(degree, perm));

    std::stringstream ss;
    ss << "Group generated correctly, generators are: ";
    for (auto j = 0u; j < generators.size(); ++j)
      ss << generators[j] << (j == generators.size() - 1u ? "" : ", ");

    EXPECT_TRUE(perm_group_equal(
      expected_elements[i], PermGroup(degree, generators, schreier_var)))
      << ss.str();
  }
}

INSTANTIATE_TEST_CASE_P(SchreierSimsVariants, SchreierSimsVariantTest,
  testing::Values(SchreierSims::SIMPLE, SchreierSims::RANDOM));
