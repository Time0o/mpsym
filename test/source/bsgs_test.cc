#include <memory>
#include <vector>

#include "gmock/gmock.h"

#include "bsgs.h"
#include "perm.h"
#include "perm_group.h"
#include "schreier_structure.h"

#include "test_main.cc"

using cgtl::BSGS;
using cgtl::Perm;
using cgtl::PermGroup;
using cgtl::SchreierTree;

TEST(BSGSTest, CanRemoveRedundantGenerators)
{
  // explicitly construct BSGS
  BSGS bsgs;

  bsgs.base = {3, 1, 2};

  bsgs.strong_generators = {
    Perm(4, {{1, 3, 4}}),
    Perm(4, {{1, 2, 3}}),
    Perm(4, {{1, 2}}),
    Perm(4, {{1, 3, 2, 4}}),
    Perm(4, {{1, 3, 2}}),
    Perm(4, {{1, 4, 3, 2}}),
    Perm(4, {{2, 3, 4}}),
    Perm(4, {{2, 3}}),
    Perm(4, {{2, 4, 3}}),
    Perm(4, {{2, 4}}),
    Perm(4, {{3, 4}})
  };

  std::vector<Perm> S1(bsgs.strong_generators);
  std::vector<Perm> S2 {Perm(4, {{1, 2}}), Perm(4, {{2, 4}})};
  std::vector<Perm> S3 {Perm(4, {{2, 4}})};

  for (int i = 0; i < 3; ++i)
    bsgs.schreier_structures.push_back(std::make_shared<SchreierTree>(4));

  bsgs.update_schreier_structure(0, S1);
  bsgs.update_schreier_structure(1, S2);
  bsgs.update_schreier_structure(2, S3);

  // remove redundant generators
  bsgs.remove_generators();

  // check whether resulting BSGS is still equivalent
  Perm expected_elements[] = {
    Perm(4, {{1, 2}}),
    Perm(4, {{1, 2}, {3, 4}}),
    Perm(4, {{1, 2, 3}}),
    Perm(4, {{1, 2, 3, 4}}),
    Perm(4, {{1, 2, 4}}),
    Perm(4, {{1, 2, 4, 3}}),
    Perm(4, {{1, 3}}),
    Perm(4, {{1, 3}, {2, 4}}),
    Perm(4, {{1, 3, 2}}),
    Perm(4, {{1, 3, 2, 4}}),
    Perm(4, {{1, 3, 4}}),
    Perm(4, {{1, 3, 4, 2}}),
    Perm(4, {{1, 4}}),
    Perm(4, {{1, 4}, {2, 3}}),
    Perm(4, {{1, 4, 2}}),
    Perm(4, {{1, 4, 2, 3}}),
    Perm(4, {{1, 4, 3}}),
    Perm(4, {{1, 4, 3, 2}}),
    Perm(4, {{2, 3}}),
    Perm(4, {{2, 3, 4}}),
    Perm(4, {{2, 4}}),
    Perm(4, {{2, 4, 3}}),
    Perm(4, {{3, 4}})
  };

  for (Perm const &perm : expected_elements) {
    EXPECT_TRUE(bsgs.strips_completely(perm))
      << "BSGS with reduced strong generators describes the same permutation group "
      << "(containing element " << perm << ")";
  }
}

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
    EXPECT_TRUE(bsgs.strips_completely(perm))
      << "Solvable group BSGS correct.";
  }

  BSGS dummy(BSGS::solve(generators_non_solvable));
  ASSERT_TRUE(dummy.base.empty())
    << "BSGS::solve fails for non-solvable group generating set.";
}
