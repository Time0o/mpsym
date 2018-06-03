#include <sstream>
#include <vector>

#include "gmock/gmock.h"
#include "perm.h"

#include "main.cc"

using testing::ElementsAre;
using testing::UnorderedElementsAre;

static ::testing::AssertionResult perm_equal(
  std::vector<unsigned> const &expected, cgtl::Perm const &perm)
{
  bool success = true;

  if (perm.degree() != expected.size()) {
    return ::testing::AssertionFailure()
      << "Permutation has incorrect degree (expected " << expected.size()
      << " but got " << perm.degree();
  }

  std::stringstream err;
  err << "Permutation differs:\n";

  for (unsigned i = 0u; i < perm.degree(); ++i) {
    if (perm[i + 1u] != expected[i]) {
      success = false;
      err << "@ index " << i + 1u << ":"
          << " expected " << expected[i]
          << " but got " << perm[i + 1u] << '\n';
    }
  }

  if (!success)
    return ::testing::AssertionFailure() << err.str();
  else
    return ::testing::AssertionSuccess();
}

TEST(PermTest, CanConstructPerm)
{
  cgtl::Perm perm;
  EXPECT_TRUE(perm_equal({1}, perm))
    << "Default construction produces identity permutation.";

  cgtl::Perm perm_id(5);
  EXPECT_TRUE(perm_equal({1, 2, 3, 4, 5}, perm_id))
    << "Identity construction produces identity permutation.";

  cgtl::Perm perm_explicit({1, 3, 4, 5, 2});
  EXPECT_TRUE(perm_equal({1, 3, 4, 5, 2}, perm_explicit))
    << "Explicit construction produces correct permutation.";

  cgtl::Perm perm_empty_cycle(6, {});
  EXPECT_TRUE(perm_equal({1, 2, 3, 4, 5, 6}, perm_empty_cycle))
    << "No-cycles construction produces correct permutation.";

  cgtl::Perm perm_single_cycle(6, {{3, 2, 5}});
  EXPECT_TRUE(perm_equal({1, 5, 2, 4, 3, 6}, perm_single_cycle))
    << "Single-cycle construction produces correct permutation.";

  cgtl::Perm perm_multi_cycles(6, {{6, 2, 4}, {2, 5, 4}, {3, 2, 5}});
  EXPECT_TRUE(perm_equal({1, 5, 2, 6, 4, 3}, perm_multi_cycles))
    << "Multi-cycle construction produces correct permutation.";
}

TEST(PermTest, CanInvertPerm)
{
  cgtl::Perm perm({3, 2, 4, 1});

  EXPECT_TRUE(perm_equal({4, 2, 1, 3}, ~perm))
    << "Inverting permutation works.";
}

TEST(PermTest, CanMultiplyPerms)
{
  cgtl::Perm perm0(7, {{1, 2, 4}});
  perm0 *= cgtl::Perm(7, {{4, 5}});

  EXPECT_TRUE(perm_equal({2, 5, 3, 1, 4, 6, 7}, perm0))
    << "Multiplying plus assigning permutation produces correct result.";

  cgtl::Perm perm1(6, {{2, 5, 4}});
  cgtl::Perm perm2(6, {{3, 2, 5}});

  cgtl::Perm perm_mult1 = perm1 * perm2;
  EXPECT_TRUE(perm_equal({1, 3, 2, 5, 4, 6}, perm_mult1))
    << "Multiplying permutations produces correct result.";
}

TEST(PermTest, PermStringRepresentation)
{
  cgtl::Perm perm1({2, 3, 1, 5, 4});
  std::stringstream ss1;
  ss1 << perm1;

  EXPECT_EQ("(1 2 3)(4 5)", ss1.str())
    << "Correct permutation string representation.";

  cgtl::Perm perm2({1, 5, 3, 6, 2, 7, 4, 8});
  std::stringstream ss2;
  ss2 << perm2;

  EXPECT_EQ("(2 5)(4 6 7)", ss2.str())
    << "Permutation string representation ignores single element cycles.";

  cgtl::Perm perm3({1, 2, 3});
  std::stringstream ss3;
  ss3 << perm3;

  EXPECT_EQ("()", ss3.str())
    << "Identity permutation string representation correct.";
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

TEST(PermGroupTest, CanTestMembership)
{
  // alternating group A4 (order 12)
  cgtl::Perm gen1(4, {{1, 2, 3}});
  cgtl::Perm gen2(4, {{1, 2}, {3, 4}});
  cgtl::PermGroup a4(4, {gen1, gen2});

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
    EXPECT_TRUE(a4.is_member(perm))
      << "Membership test correctly identifies group member " << perm;
  }

  for (cgtl::Perm const &perm : expected_non_members) {
    EXPECT_FALSE(a4.is_member(perm))
      << "Membership test correctly rejects non group member " << perm;
  }
}
