#include <algorithm>
#include <string>
#include <sstream>

#include "gmock/gmock.h"

#include "partial_perm.h"

#include "test_main.cc"

using cgtl::PartialPerm;

TEST(PartialPermTest, CanConstructPartialPerm)
{
  std::vector<std::vector<unsigned>> pperms {
    {0, 4, 0, 3, 0, 9, 6, 0, 7, 0, 11},
    {5, 9, 10, 11, 0, 0, 0, 0, 0, 12, 4, 3}
  };

  std::vector<std::vector<unsigned>> expected_doms {
    {2, 4, 6, 7, 9, 11},
    {1, 2, 3, 4, 10, 11, 12}
  };

  std::vector<std::vector<unsigned>> expected_ims {
    {3, 4, 6, 7, 9, 11},
    {3, 4, 5, 9, 10, 11, 12}
  };

  for (auto i = 0u; i < pperms.size(); ++i) {
    PartialPerm partial_perm(pperms[i]);

    for (auto j = 0u; j < pperms[i].size(); ++j) {
      EXPECT_EQ(pperms[i][j], partial_perm[j + 1u])
        << "Can apply partial permutation.";
    }

    EXPECT_EQ(expected_doms[i], partial_perm.dom())
      << "Partial permutation domain constructed correct.";

    EXPECT_EQ(*std::min_element(expected_doms[i].begin(), expected_doms[i].end()),
              partial_perm.dom_min())
      << "Partial permutation domain lower limit correct.";

    EXPECT_EQ(*std::max_element(expected_doms[i].begin(), expected_doms[i].end()),
              partial_perm.dom_max())
      << "Partial permutation domain upper limit correct.";

    EXPECT_EQ(expected_ims[i], partial_perm.im())
      << "Partial permutation image constructed correct.";

    EXPECT_EQ(*std::min_element(expected_ims[i].begin(), expected_ims[i].end()),
              partial_perm.im_min())
      << "Partial permutation image lower limit correct.";

    EXPECT_EQ(*std::max_element(expected_ims[i].begin(), expected_ims[i].end()),
              partial_perm.im_max())
      << "Partial permutation image upper limit correct.";
  }
}

TEST(PartialPermTest, CanInvertPartialPerm)
{
  PartialPerm inv(~PartialPerm({0, 4, 0, 3, 0, 9, 6, 0, 7, 0, 11}));
  PartialPerm expected({0, 0, 4, 2, 0, 7, 9, 0, 6, 0, 11});

  ASSERT_EQ(expected, inv)
    << "Inverting partial permutation produces correct result.";

  ASSERT_EQ(expected.dom(), inv.dom())
    << "Inverting partial permutation produces correct domain.";

  ASSERT_EQ(expected.im(), inv.im())
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
    ASSERT_EQ(expected, pperm)
      << "Multiplying partial permutations produces correct result.";

    ASSERT_EQ(expected.dom(), pperm.dom())
      << "Multiplying partial permutations produces correct domain.";

    ASSERT_EQ(expected.im(), pperm.im())
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
