#include <sstream>
#include <vector>

#include "perm.h"
#include "gtest/gtest.h"

enum { DEFAULT_N = 10u };

template<unsigned N>
static ::testing::AssertionResult perm_equal(
  std::vector<unsigned> const &expected, CGTL::Perm<N> const &perm)
{
  bool success = true;

  std::stringstream err;
  err << "Permutation differs:\n";

  for (unsigned i = 0u; i < N; ++i) {
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
  CGTL::Perm<DEFAULT_N> perm_id;
  for (unsigned i = 1u; i <= DEFAULT_N; ++i) {
    EXPECT_EQ(i, perm_id[i])
      << "Default constructor produces identity permutation.";
  }

  CGTL::Perm<6> perm_cycle({ 3, 2, 5 });
  EXPECT_TRUE(perm_equal({ 1, 5, 2, 4, 3, 6 }, perm_cycle))
    << "Single-cycle constructor produces correct permutation.";

  CGTL::Perm<6> perm_multicycles({ { 6, 2, 4}, { 2, 5, 4 }, { 3, 2, 5 } });
  EXPECT_TRUE(perm_equal({ 1, 6, 5, 4, 3, 2 }, perm_multicycles))
    << "Multi-cycle constructor produces correct permutation.";

  CGTL::Perm<10> perm_widened1(perm_cycle);
  CGTL::Perm<10> perm_widened2 = perm_cycle;

  EXPECT_TRUE(perm_equal({ 1, 5, 2, 4, 3, 6, 7, 8, 9, 10 }, perm_widened1))
     << "Widening permutation works (direct initialization).";
  EXPECT_TRUE(perm_equal({ 1, 5, 2, 4, 3, 6, 7, 8, 9, 10 }, perm_widened2))
     << "Widening permutation works (copy initialization).";
}

TEST(PermTest, CanMultiplyPerms)
{
  CGTL::Perm<6> perm1({ 2, 5, 4 });
  CGTL::Perm<6> perm2({ 3, 2, 5 });

  CGTL::Perm<6> perm_mult1 = perm1 * perm2;
  EXPECT_TRUE(perm_equal({ 1, 4, 5, 2, 3, 6 }, perm_mult1))
    << "Multiplying permutations produces correct result.";

  CGTL::Perm<3> perm3({ 1, 2 });
  CGTL::Perm<6> perm4({ 4, 1, 5, 2 });

  CGTL::Perm<6> perm_mult2 = perm3 * perm4;
  EXPECT_TRUE(perm_equal({ 5, 4, 3, 2, 1, 6 }, perm_mult2))
     << "Right multiplying larger permutation produces correct result.";

  CGTL::Perm<6> perm_mult3 = perm4 * perm3;
  EXPECT_TRUE(perm_equal({ 4, 5, 3, 1, 2, 6 }, perm_mult3))
     << "Left multiplying larger permutation produces correct result.";
}
