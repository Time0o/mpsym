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
  EXPECT_TRUE(perm_equal({1}, perm))
    << "Default construction produces identity permutation.";

  Perm perm_id(5);
  EXPECT_TRUE(perm_equal({1, 2, 3, 4, 5}, perm_id))
    << "Identity construction produces identity permutation.";

  Perm perm_explicit({1, 3, 4, 5, 2});
  EXPECT_TRUE(perm_equal({1, 3, 4, 5, 2}, perm_explicit))
    << "Explicit construction produces correct permutation.";

  Perm perm_empty_cycle(6, {});
  EXPECT_TRUE(perm_equal({1, 2, 3, 4, 5, 6}, perm_empty_cycle))
    << "No-cycles construction produces correct permutation.";

  Perm perm_single_cycle(6, {{3, 2, 5}});
  EXPECT_TRUE(perm_equal({1, 5, 2, 4, 3, 6}, perm_single_cycle))
    << "Single-cycle construction produces correct permutation.";

  Perm perm_multi_cycles(6, {{6, 2, 4}, {2, 5, 4}, {3, 2, 5}});
  EXPECT_TRUE(perm_equal({1, 5, 2, 6, 4, 3}, perm_multi_cycles))
    << "Multi-cycle construction produces correct permutation.";
}

TEST(PermTest, CanInvertPerm)
{
  Perm perm({3, 2, 4, 1});

  EXPECT_TRUE(perm_equal({4, 2, 1, 3}, ~perm))
    << "Inverting permutation works.";
}

TEST(PermTest, CanMultiplyPerms)
{
  Perm perm0(7, {{1, 2, 4}});
  perm0 *= Perm(7, {{4, 5}});

  EXPECT_TRUE(perm_equal({2, 5, 3, 1, 4, 6, 7}, perm0))
    << "Multiplying plus assigning permutation produces correct result.";

  Perm perm1(6, {{2, 5, 4}});
  Perm perm2(6, {{3, 2, 5}});

  Perm perm_mult1 = perm1 * perm2;
  EXPECT_TRUE(perm_equal({1, 3, 2, 5, 4, 6}, perm_mult1))
    << "Multiplying permutations produces correct result.";
}

TEST(PermTest, PermStringRepresentation)
{
  Perm perm1({2, 3, 1, 5, 4});
  std::stringstream ss1;
  ss1 << perm1;

  EXPECT_EQ("(1, 2, 3)(4, 5)", ss1.str())
    << "Correct permutation string representation.";

  Perm perm2({1, 5, 3, 6, 2, 7, 4, 8});
  std::stringstream ss2;
  ss2 << perm2;

  EXPECT_EQ("(2, 5)(4, 6, 7)", ss2.str())
    << "Permutation string representation ignores single element cycles.";

  Perm perm3({1, 2, 3});
  std::stringstream ss3;
  ss3 << perm3;

  EXPECT_EQ("()", ss3.str())
    << "Identity permutation string representation correct.";
}

TEST(PermTest, CanHashPerm)
{
  std::vector<Perm> perms = {
    Perm(5, {{1, 2, 3}}),
    Perm(5, {{2, 3}, {4, 5}}),
    Perm(5, {{1, 2, 3, 4}}),
    Perm(5, {{1, 2}}),
    Perm(5, {{1, 2, 3}, {4, 5}})
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
  Perm perm(5, {{2, 5}, {3, 1, 4}});

  std::vector<std::vector<unsigned>> expected_extended_perms {
    {4, 5, 1, 3, 2},
    {4, 5, 1, 3, 2, 6},
    {4, 5, 1, 3, 2, 6, 7},
    {4, 5, 1, 3, 2, 6, 7, 8}
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
  Perm perm(5, {{2, 5}, {3, 1, 4}});

  std::vector<std::vector<unsigned>> expected_shifted_perms {
    {4, 5, 1, 3, 2},
    {1, 5, 6, 2, 4, 3},
    {1, 2, 6, 7, 3, 5, 4},
    {1, 2, 3, 7, 8, 4, 6, 5},
    {1, 2, 3, 4, 8, 9, 5, 7, 6},
    {1, 2, 3, 4, 5, 9, 10, 6, 8, 7},
    {1, 2, 3, 4, 5, 6, 10, 11, 7, 9, 8},
    {1, 2, 3, 4, 5, 6, 7, 11, 12, 8, 10, 9}
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
      Perm(4, {{1, 2, 3}}),
      {1, 2, 3},
      Perm(4, {{1, 2, 3}})),
    PermRestriction(
      Perm(9, {{2, 4}, {3, 5}, {1, 7, 8}}),
      {2, 3, 4, 5},
      Perm(9, {{2, 4}, {3, 5}})),
    PermRestriction(
      Perm(12, {{6, 3, 2, 1}, {4, 7}, {9, 8}, {10, 11}}),
      {4, 7, 8, 9, 10, 11},
      Perm(12, {{4, 7}, {8, 9}, {10, 11}})),
    PermRestriction(
      Perm(3, {{1, 2}}),
      {3},
      Perm(3)),
    PermRestriction(
      Perm(5, {{1, 2, 3}}),
      {4, 5},
      Perm(5))
  };

  for (auto const &pr : perm_restrictions) {
    EXPECT_EQ(pr.expected, pr.perm.restricted(pr.domain.begin(), pr.domain.end()))
      << "Restricting permutation yields correct result.";
  }
}
