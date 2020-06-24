#include <memory>
#include <vector>

#include "gmock/gmock.h"

#include "bsgs.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"
#include "schreier_structure.h"

#include "test_main.cc"

using mpsym::BSGS;
using mpsym::BSGSOptions;
using mpsym::Perm;
using mpsym::PermGroup;
using mpsym::PermSet;

TEST(DISABLED_BSGSSolveTest, CanSolveBSGS)
{
  BSGSOptions bsgs_options;
  bsgs_options.construction = BSGSOptions::Construction::SOLVE;

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

  PermSet generators_non_solvable(PermGroup::symmetric(5).generators());

  BSGS bsgs(4, generators_solvable, &bsgs_options);

  for (Perm const &perm : generators_solvable_expected_elements) {
    EXPECT_TRUE(bsgs.strips_completely(perm))
      << "Solvable group BSGS correct.";
  }

  EXPECT_THROW(BSGS dummy(5, generators_non_solvable, &bsgs_options),
               BSGS::SolveError)
      << "Solving BSGS fails for non-solvable group generating set.";
}

//TEST(BSGSBaseSwapTest, CanConjugateBSGS)
//{
//  PermGroup pg(5, {Perm(5, {{1, 2}, {3, 4}}), Perm(5, {{1, 4, 2}})});
//
//  for (Perm const &perm : pg) {
//    PermGroup pg_conjugated_bsgs(pg);
//    pg_conjugated_bsgs.bsgs().conjugate(perm)
//
//    EXPECT_TRUE(perm_group_equal(pg, pg_conjugated_bsgs))
//      << "Permutation group remains the same after conjugating BSGS."
//  }
//}
