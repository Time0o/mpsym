#include <algorithm>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "gmock/gmock.h"

#include "partial_perm.h"
#include "perm.h"
#include "test_utility.h"

#include "test_main.cc"

using cgtl::PartialPerm;
using cgtl::Perm;

TEST(PartialPermTest, CanConstructPartialPerm)
{
  std::vector<PartialPerm> pperms {
    PartialPerm(),
    PartialPerm({}),
    PartialPerm({}, {}),
    PartialPerm::id({}),
    PartialPerm(5),
    PartialPerm::id({3, 5, 4}),
    PartialPerm({0, 4, 0, 3, 0, 9, 6, 0, 7, 0, 11}),
    PartialPerm({2, 4, 6, 7, 9, 11}, {4, 3, 9, 6, 7, 11}),
    PartialPerm({5, 9, 10, 11, 0, 0, 0, 0, 0, 12, 4, 3}),
    PartialPerm({12, 11, 1, 2, 3, 4, 10}, {3, 4, 5, 9, 10, 11, 12})
  };

  std::vector<std::vector<unsigned>> expected_mappings {
    {},
    {},
    {},
    {},
    {1, 2, 3, 4, 5},
    {0, 0, 3, 4, 5},
    {0, 4, 0, 3, 0, 9, 6, 0, 7, 0, 11},
    {0, 4, 0, 3, 0, 9, 6, 0, 7, 0, 11},
    {5, 9, 10, 11, 0, 0, 0, 0, 0, 12, 4, 3},
    {5, 9, 10, 11, 0, 0, 0, 0, 0, 12, 4, 3}
  };

  std::vector<std::vector<unsigned>> expected_doms {
    {},
    {},
    {},
    {},
    {1, 2, 3, 4, 5},
    {3, 4, 5},
    {2, 4, 6, 7, 9, 11},
    {2, 4, 6, 7, 9, 11},
    {1, 2, 3, 4, 10, 11, 12},
    {1, 2, 3, 4, 10, 11, 12}
  };

  std::vector<std::vector<unsigned>> expected_ims {
    {},
    {},
    {},
    {},
    {1, 2, 3, 4, 5},
    {3, 4, 5},
    {3, 4, 6, 7, 9, 11},
    {3, 4, 6, 7, 9, 11},
    {3, 4, 5, 9, 10, 11, 12},
    {3, 4, 5, 9, 10, 11, 12}
  };

  for (auto i = 0u; i < pperms.size(); ++i) {
    for (auto j = 0u; j < expected_mappings[i].size(); ++j) {
      EXPECT_EQ(expected_mappings[i][j], pperms[i][j + 1u])
        << "Can apply partial permutation.";
    }

    EXPECT_EQ(expected_doms[i], pperms[i].dom())
      << "Partial permutation domain constructed correct.";

    if (expected_doms[i].empty()) {
      EXPECT_EQ(0u, pperms[i].dom_min())
        << "Partial permutation domain lower limit correct.";

      EXPECT_EQ(0u, pperms[i].dom_max())
        << "Partial permutation domain uppter limit correct.";

    } else {
      EXPECT_EQ(*std::min_element(expected_doms[i].begin(), expected_doms[i].end()),
                pperms[i].dom_min())
        << "Partial permutation domain lower limit correct.";

      EXPECT_EQ(*std::max_element(expected_doms[i].begin(), expected_doms[i].end()),
                pperms[i].dom_max())
        << "Partial permutation domain upper limit correct.";
    }

    EXPECT_EQ(expected_ims[i], pperms[i].im())
      << "Partial permutation image constructed correct.";

    if (expected_ims[i].empty()) {
      EXPECT_EQ(0u, pperms[i].im_min())
        << "Partial permutation domain lower limit correct.";

      EXPECT_EQ(0u, pperms[i].im_max())
        << "Partial permutation domain uppter limit correct.";

    } else {
      EXPECT_EQ(*std::min_element(expected_ims[i].begin(), expected_ims[i].end()),
                pperms[i].im_min())
        << "Partial permutation image lower limit correct.";

      EXPECT_EQ(*std::max_element(expected_ims[i].begin(), expected_ims[i].end()),
                pperms[i].im_max())
        << "Partial permutation image upper limit correct.";
    }
  }
}

TEST(PartialPermTest, CanInvertPartialPerm)
{
  PartialPerm inv(~PartialPerm({0, 4, 0, 3, 0, 9, 6, 0, 7, 0, 11}));
  PartialPerm expected({0, 0, 4, 2, 0, 7, 9, 0, 6, 0, 11});

  EXPECT_EQ(expected, inv)
    << "Inverting partial permutation produces correct result.";

  EXPECT_EQ(expected.dom(), inv.dom())
    << "Inverting partial permutation produces correct domain.";

  EXPECT_EQ(expected.im(), inv.im())
    << "Inverting partial permutation produces correct image.";

  EXPECT_TRUE(expected.dom_min() == inv.dom_min() &&
              expected.dom_max() == inv.dom_max())
    << "Inverting partial permutation produces correct domain limits.";

  EXPECT_TRUE(expected.im_min() == inv.im_min() &&
              expected.im_max() == inv.im_max())
    << "Inverting partial permutation produces correct image limits.";
}

TEST(PartialPermTest, CanMultiplyPartialPerms)
{
  PartialPerm lhs({0, 4, 0, 3, 0, 9, 6, 0, 7, 0, 11});
  PartialPerm rhs({5, 9, 10, 11, 0, 0, 0, 0, 0, 12, 4, 3});
  PartialPerm expected({0, 11, 0, 10, 0, 0, 0, 0, 0, 0, 4});

  PartialPerm pperm_mult_assign(lhs);
  pperm_mult_assign *= rhs;

  PartialPerm pperm_mult = lhs * rhs;

  for (PartialPerm const &pperm : {pperm_mult_assign, pperm_mult}) {
    EXPECT_EQ(expected, pperm)
      << "Multiplying partial permutations produces correct result.";

    EXPECT_EQ(expected.dom(), pperm.dom())
      << "Multiplying partial permutations produces correct domain.";

    EXPECT_EQ(expected.im(), pperm.im())
      << "Multiplying partial permutations produces correct image.";

    EXPECT_TRUE(expected.dom_min() == pperm.dom_min() &&
                expected.dom_max() == pperm.dom_max())
      << "Multiplying partial permutations produces correct domain limits.";

    EXPECT_TRUE(expected.im_min() == pperm.im_min() &&
                expected.im_max() == pperm.im_max())
      << "Multiplying partial permutations produces correct image limits.";
  }
}

TEST(PartialPermTest, PartialPermStringRepresentation)
{
  struct PPermStrRepr {
    PPermStrRepr(PartialPerm const &pperm, std::string const &str)
      : pperm(pperm), str(str) {}

    PartialPerm pperm;
    std::string str;
  };

  PPermStrRepr pperm_str_reprs[] = {
    PPermStrRepr(PartialPerm(),
                 "()"),
    PPermStrRepr(PartialPerm({1, 0, 3}),
                 "(1)(3)"),
    PPermStrRepr(PartialPerm({0, 2, 0}),
                 "(2)"),
    PPermStrRepr(PartialPerm({2, 0, 0, 1}),
                 "[4 1 2]"),
    PPermStrRepr(PartialPerm({0, 1, 5, 0, 2}),
                 "[3 5 2 1]"),
    PPermStrRepr(PartialPerm({0, 0, 3, 4, 1, 0}),
                 "[5 1](3)(4)"),
    PPermStrRepr(PartialPerm({6, 9, 7, 1, 0, 5, 3, 10, 0, 11, 8}),
                 "[2 9][4 1 6 5](3 7)(8 10 11)")
  };

  for (auto const &pperm_str_repr : pperm_str_reprs) {
    std::stringstream ss;
    ss << pperm_str_repr.pperm;

    EXPECT_EQ(pperm_str_repr.str, ss.str())
      << "Correct partial permutation string representation.";
  }
}

TEST(PartialPermTest, CanRestrictPartialPerm)
{
  struct RestrictionTest {
    RestrictionTest(PartialPerm pperm, std::vector<unsigned> const &domain,
                    PartialPerm expected)
      : pperm(pperm), domain(domain), expected(expected) {}

    PartialPerm pperm;
    std::vector<unsigned> domain;
    PartialPerm expected;
  };

  RestrictionTest tests[] = {
    RestrictionTest(
      PartialPerm({0, 4, 0, 3, 0, 9, 6, 0, 7, 0, 11}),
      {4, 5, 6, 9, 10},
      PartialPerm({0, 0, 0, 3, 0, 9, 0, 0, 7})),
    RestrictionTest(
      PartialPerm({5, 9, 10, 11, 0, 0, 0, 0, 0, 12, 4, 3}),
      {1, 2, 3, 8, 9},
      PartialPerm({5, 9, 10}))
  };

  for (auto &test : tests) {
    PartialPerm actual(test.pperm.restricted(test.domain));

    EXPECT_EQ(test.expected, actual)
      << "Multiplying partial permutations produces correct result.";

    EXPECT_EQ(test.expected.dom(), actual.dom())
      << "Multiplying partial permutations produces correct domain.";

    EXPECT_EQ(test.expected.im(), actual.im())
      << "Multiplying partial permutations produces correct image.";

    EXPECT_TRUE(test.expected.dom_min() == actual.dom_min() &&
                test.expected.dom_max() == actual.dom_max())
      << "Multiplying partial permutations produces correct domain limits.";

    EXPECT_TRUE(test.expected.im_min() == actual.im_min() &&
                test.expected.im_max() == actual.im_max())
      << "Multiplying partial permutations produces correct image limits.";
  }
}

TEST(PartialPermTest, CanConvertPartialPermToPerm)
{
  std::vector<std::pair<PartialPerm, Perm>> conversions {
    {
      PartialPerm(),
      Perm()
    },
    {
      PartialPerm(),
      Perm(10, {})
    },
    {
      PartialPerm({1, 2}, {2, 1}),
      Perm(3, {{1, 2}})
    },
    {
      PartialPerm({2, 3, 5}, {3, 2, 5}),
      Perm(6, {{2, 3}})
    },
    {
      PartialPerm({4, 5, 6, 7, 8, 9}, {4, 7, 8, 5, 9, 6}),
      Perm(10, {{5, 7}, {6, 8, 9}})
    }
  };

  for (auto const &conv : conversions) {
    PartialPerm const &pperm = std::get<0>(conv);
    Perm const &perm = std::get<1>(conv);

    EXPECT_EQ(perm, pperm.to_perm(perm.degree()))
      << "Conversion from partial to 'complete' permutation correct.";
  }
}
