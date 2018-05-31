#include <sstream>
#include <vector>

#include "perm.h"
#include "gtest/gtest.h"

static ::testing::AssertionResult perm_equal(
  std::vector<unsigned> const &expected, cgtl::Perm const &perm)
{
  bool success = true;

  std::stringstream err;
  err << "Permutation differs:\n";

  for (unsigned i = 0u; i < perm.size(); ++i) {
    if (perm[i + 1] != expected[i]) {
      success = false;
      err << "@ index " << i + 1 << ":"
          << " expected " << expected[i]
          << " but got " << perm[i + 1] << '\n';
    }
  }

  if (!success)
    return ::testing::AssertionFailure() << err.str();
  else
    return ::testing::AssertionSuccess();
}

TEST(PermTest, CanConstructPerm)
{
  cgtl::Perm perm_id(10);
  for (unsigned i = 1; i <= 10; ++i) {
    EXPECT_EQ(i, perm_id[i])
      << "Default constructor produces identity permutation.";
  }

  cgtl::Perm perm_explicit({ 1, 3, 4, 5, 2 });
  EXPECT_EQ(5u, perm_explicit.size())
    << "Explicit construction produces permutation of correct size.";
  EXPECT_TRUE(perm_equal({ 1, 3, 4, 5, 2 }, perm_explicit))
    << "Explicit construction produces correct permutation.";

  cgtl::Perm perm_empty_cycle(6, { });
  EXPECT_TRUE(perm_equal({ 1, 2, 3, 4, 5, 6 }, perm_empty_cycle))
    << "No-cycles construction produces correct permutation.";

  cgtl::Perm perm_single_cycle(6, { { 3, 2, 5 } });
  EXPECT_TRUE(perm_equal({ 1, 5, 2, 4, 3, 6 }, perm_single_cycle))
    << "Single-cycle construction produces correct permutation.";

  cgtl::Perm perm_multi_cycles(6, { { 6, 2, 4 }, { 2, 5, 4 }, { 3, 2, 5 } });
  EXPECT_TRUE(perm_equal({ 1, 6, 5, 4, 3, 2 }, perm_multi_cycles))
    << "Multi-cycle construction produces correct permutation.";
}

TEST(PermTest, CanExtendPerm)
{
  cgtl::Perm perm{3, { { 1, 2 } }};

  EXPECT_TRUE(perm_equal({ 2, 1, 3, 4, 5, 6 }, perm.extend(6)))
     << "Extending permutation works.";
}

TEST(PermTest, CanMultiplyPerms)
{
  cgtl::Perm perm1(6, { { 2, 5, 4 } });
  cgtl::Perm perm2(6, { { 3, 2, 5 } });

  cgtl::Perm perm_mult1 = perm1 * perm2;
  EXPECT_TRUE(perm_equal({ 1, 4, 5, 2, 3, 6 }, perm_mult1))
    << "Multiplying permutations produces correct result.";

  cgtl::Perm perm3(3, { { 1, 2 } });
  cgtl::Perm perm4(6, { { 4, 1, 5, 2 } });

  cgtl::Perm perm_mult2 = perm3 * perm4;
  EXPECT_EQ(6u, perm_mult2.size())
     << "Right multiplying larger permutation produces result of correct size.";
  EXPECT_TRUE(perm_equal({ 5, 4, 3, 2, 1, 6 }, perm_mult2))
     << "Right multiplying larger permutation produces correct result.";

  cgtl::Perm perm_mult3 = perm4 * perm3;
  EXPECT_EQ(6u, perm_mult3.size())
     << "Left multiplying larger permutation produces result of correct size.";
  EXPECT_TRUE(perm_equal({ 4, 5, 3, 1, 2, 6 }, perm_mult3))
     << "Left multiplying larger permutation produces correct result.";
}
