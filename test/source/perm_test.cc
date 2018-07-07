#include <vector>

#include "gmock/gmock.h"
#include "perm.h"
#include "test_utility.h"

#include "test_main.cc"

using cgtl::Perm;

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

  EXPECT_EQ("(1 2 3)(4 5)", ss1.str())
    << "Correct permutation string representation.";

  Perm perm2({1, 5, 3, 6, 2, 7, 4, 8});
  std::stringstream ss2;
  ss2 << perm2;

  EXPECT_EQ("(2 5)(4 6 7)", ss2.str())
    << "Permutation string representation ignores single element cycles.";

  Perm perm3({1, 2, 3});
  std::stringstream ss3;
  ss3 << perm3;

  EXPECT_EQ("()", ss3.str())
    << "Identity permutation string representation correct.";
}
