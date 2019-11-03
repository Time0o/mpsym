#include <memory>
#include <vector>

#include "gmock/gmock.h"

#include "bsgs.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"
#include "schreier_structure.h"

#include "test_main.cc"

using cgtl::BSGS;
using cgtl::Perm;
using cgtl::PermGroup;
using cgtl::PermSet;
using cgtl::SchreierTree;

TEST(BSGSTest, CanSolveBSGS)
{
  PermSet generators_solvable {
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

  PermSet generators_non_solvable(
    PermGroup::symmetric(5).bsgs().strong_generators());

  BSGS bsgs(4, generators_solvable, BSGS::CONSTRUCTION_SOLVE);

  for (Perm const &perm : generators_solvable_expected_elements) {
    EXPECT_TRUE(bsgs.strips_completely(perm))
      << "Solvable group BSGS correct.";
  }

  EXPECT_THROW(BSGS dummy(5, generators_non_solvable, BSGS::CONSTRUCTION_SOLVE),
               BSGS::SolveError)
      << "Solving BSGS fails for non-solvable group generating set.";
}
