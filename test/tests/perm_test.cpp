#include <sstream>
#include <unordered_set>
#include <vector>

#include "gmock/gmock.h"

#include "perm.hpp"
#include "test_utility.hpp"

#include "test_main.cpp"

using namespace mpsym;
using namespace mpsym::internal;

using testing::UnorderedElementsAreArray;

TEST(PermTest, CanConstructPerm)
{
  Perm perm;
  EXPECT_TRUE(perm_equal({0}, perm))
    << "Default construction produces identity permutation.";

  Perm perm_id(5);
  EXPECT_TRUE(perm_equal({0, 1, 2, 3, 4}, perm_id))
    << "Identity construction produces identity permutation.";

  Perm perm_explicit({0, 2, 3, 4, 1});
  EXPECT_TRUE(perm_equal({0, 2, 3, 4, 1}, perm_explicit))
    << "Explicit construction produces correct permutation.";

  Perm perm_empty_cycle(6, {});
  EXPECT_TRUE(perm_equal({0, 1, 2, 3, 4, 5}, perm_empty_cycle))
    << "No-cycles construction produces correct permutation.";

  Perm perm_single_cycle(6, {{2, 1, 4}});
  EXPECT_TRUE(perm_equal({0, 4, 1, 3, 2, 5}, perm_single_cycle))
    << "Single-cycle construction produces correct permutation.";

  Perm perm_multi_cycles(6, {{5, 1, 3}, {1, 4, 3}, {2, 1, 4}});
  EXPECT_TRUE(perm_equal({0, 4, 1, 5, 3, 2}, perm_multi_cycles))
    << "Multi-cycle construction produces correct permutation.";
}

TEST(PermTest, CanInvertPerm)
{
  Perm perm({2, 1, 3, 0});

  EXPECT_TRUE(perm_equal({3, 1, 0, 2}, ~perm))
    << "Inverting permutation works.";
}

TEST(PermTest, CanMultiplyPerms)
{
  Perm perm0(7, {{0, 1, 3}});
  perm0 *= Perm(7, {{3, 4}});

  EXPECT_TRUE(perm_equal({1, 4, 2, 0, 3, 5, 6}, perm0))
    << "Multiplying plus assigning permutation produces correct result.";

  Perm perm1(6, {{1, 4, 3}});
  Perm perm2(6, {{2, 1, 4}});

  Perm perm_mult1 = perm1 * perm2;
  EXPECT_TRUE(perm_equal({0, 2, 1, 4, 3, 5}, perm_mult1))
    << "Multiplying permutations produces correct result.";
}

TEST(PermTest, PermStringRepresentation)
{
  Perm perm1({1, 2, 0, 4, 3});
  std::stringstream ss1;
  ss1 << perm1;

  EXPECT_EQ("(0, 1, 2)(3, 4)", ss1.str())
    << "Correct permutation string representation.";

  Perm perm2({0, 4, 2, 5, 1, 6, 3, 7});
  std::stringstream ss2;
  ss2 << perm2;

  EXPECT_EQ("(1, 4)(3, 5, 6)", ss2.str())
    << "Permutation string representation ignores single element cycles.";

  Perm perm3({0, 1, 2});
  std::stringstream ss3;
  ss3 << perm3;

  EXPECT_EQ("()", ss3.str())
    << "Identity permutation string representation correct.";
}

TEST(PermTest, CanHashPerm)
{
  std::vector<Perm> perms = {
    Perm(5, {{0, 1, 2}}),
    Perm(5, {{1, 2}, {3, 4}}),
    Perm(5, {{0, 1, 2, 3}}),
    Perm(5, {{0, 1}}),
    Perm(5, {{0, 1, 2}, {3, 4}})
  };

  std::unordered_set<Perm> permset;

  unsigned const repetitions = 10u;
  for (unsigned i = 0u; i < repetitions; ++i) {
    for (Perm const &p : perms)
      permset.insert(p);
  }

  ASSERT_EQ(perms.size(), permset.size())
    << "Hashed permutation set has correct size.";

  std::vector<Perm>hashed_perms(permset.begin(), permset.end());
  EXPECT_THAT(hashed_perms, UnorderedElementsAreArray(perms))
    << "Hashed permutation set has correct elements.";
}

TEST(PermTest, CanExtendPerm)
{
  Perm perm(5, {{1, 4}, {2, 0, 3}});

  std::vector<std::vector<unsigned>> expected_extended_perms {
    {3, 4, 0, 2, 1},
    {3, 4, 0, 2, 1, 5},
    {3, 4, 0, 2, 1, 5, 6},
    {3, 4, 0, 2, 1, 5, 6, 7}
  };

  for (auto i = 0u; i < expected_extended_perms.size(); ++i) {
    EXPECT_TRUE(perm_equal(expected_extended_perms[i],
                           perm.extended(perm.degree() + i)))
      << "Permutation extension yields original permutation.";
  }
}

// TODO: CanNormalizePerm

TEST(PermTest, CanShiftPerm)
{
  Perm perm(5, {{1, 4}, {2, 0, 3}});

  std::vector<std::vector<unsigned>> expected_shifted_perms {
    {3, 4, 0, 2, 1},
    {0, 4, 5, 1, 3, 2},
    {0, 1, 5, 6, 2, 4, 3},
    {0, 1, 2, 6, 7, 3, 5, 4},
    {0, 1, 2, 3, 7, 8, 4, 6, 5},
    {0, 1, 2, 3, 4, 8, 9, 5, 7, 6},
    {0, 1, 2, 3, 4, 5, 9, 10, 6, 8, 7},
    {0, 1, 2, 3, 4, 5, 6, 10, 11, 7, 9, 8}
  };

  for (auto i = 0u; i < expected_shifted_perms.size(); ++i) {
    EXPECT_TRUE(perm_equal(expected_shifted_perms[i], perm.shifted(i)))
      << "Permutation shift yields original permutation (shift was " << i << ")";
  }
}

TEST(PermTest, CanRestrictPerm)
{
  struct PermRestriction {
    PermRestriction(
      Perm const &perm, std::vector<unsigned> const &domain,
      Perm const &expected) : perm(perm), domain(domain), expected(expected) {}

    Perm perm;
    std::vector<unsigned> domain;
    Perm expected;
  };

  PermRestriction perm_restrictions[] = {
    PermRestriction(
      Perm(4, {{0, 1, 2}}),
      {0, 1, 2},
      Perm(4, {{0, 1, 2}})),
    PermRestriction(
      Perm(9, {{1, 3}, {2, 4}, {0, 6, 7}}),
      {1, 2, 3, 4},
      Perm(9, {{1, 3}, {2, 4}})),
    PermRestriction(
      Perm(12, {{5, 2, 1, 0}, {3, 6}, {8, 7}, {9, 10}}),
      {3, 6, 7, 8, 9, 10},
      Perm(12, {{3, 6}, {7, 8}, {9, 10}})),
    PermRestriction(
      Perm(3, {{0, 1}}),
      {2},
      Perm(3)),
    PermRestriction(
      Perm(5, {{0, 1, 2}}),
      {3, 4},
      Perm(5))
  };

  for (auto const &pr : perm_restrictions) {
    EXPECT_EQ(pr.expected, pr.perm.restricted(pr.domain.begin(), pr.domain.end()))
      << "Restricting permutation yields correct result.";
  }
}
