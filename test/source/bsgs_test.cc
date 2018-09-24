#include <vector>

#include "gmock/gmock.h"

#include "bsgs.h"
#include "perm.h"
#include "perm_group.h"

#include "test_main.cc"

using cgtl::BSGS;
using cgtl::Perm;
using cgtl::PermGroup;

TEST(BSGSTest, CanSolveBSGS)
{
  std::vector<Perm> generators_solvable {
    Perm(4, {{2, 4}}),
    Perm(4, {{1, 2}, {3, 4}})
  };

  Perm generators_solvable_expected_elements[] = {
    Perm(4, {{1, 2, 3, 4}}),
    Perm(4, {{1, 2}, {3, 4}}),
    Perm(4, {{1, 3}, {2, 4}}),
    Perm(4, {{1, 3}}),
    Perm(4, {{1, 4, 3, 2}}),
    Perm(4, {{1, 4}, {2, 3}}),
    Perm(4, {{2, 4}})
  };

  std::vector<Perm> generators_non_solvable(
    PermGroup::symmetric(5).bsgs().strong_generators);

  BSGS bsgs(BSGS::solve(generators_solvable));
  ASSERT_FALSE(bsgs.base.empty())
    << "BSGS::solve succeeds for solvable group generating set.";

  for (Perm const &perm : generators_solvable_expected_elements) {
    EXPECT_TRUE(bsgs.contains(perm))
      << "Solvable group BSGS correct.";
  }

  BSGS dummy(BSGS::solve(generators_non_solvable));
  ASSERT_TRUE(dummy.base.empty())
    << "BSGS::solve fails for non-solvable group generating set.";
}
