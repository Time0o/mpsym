#include <string>
#include <sstream>

#include "gmock/gmock.h"

#include "partial_perm.h"

#include "test_main.cc"

using cgtl::PartialPerm;

TEST(PartialPermTest, PartialPermStringRepresentation)
{
  struct PPermStrRepr {
    PPermStrRepr(PartialPerm const &pperm, std::string const &str)
      : pperm(pperm), str(str) {}

    PartialPerm pperm;
    std::string str;
  };

  PPermStrRepr pperm_str_reprs[] = {
    PPermStrRepr(PartialPerm({}),
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
