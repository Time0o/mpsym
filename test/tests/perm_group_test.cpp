#include <algorithm>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "gmock/gmock.h"

#include "bsgs.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"
#include "test_utility.hpp"
#include "util.hpp"

#include "test_main.cpp"

using namespace mpsym;
using namespace mpsym::internal;

using testing::ElementsAre;
using testing::UnorderedElementsAre;
using testing::UnorderedElementsAreArray;

TEST(PermGroupTest, CanCompare)
{
  PermGroup pg1(
    {
      Perm(5, {{0, 1}, {2, 3}}),
      Perm(5, {{0, 3, 1}})
    }
  );

  PermGroup pg2(
    {
      Perm(5, {{0, 1}, {2, 3}}),
      Perm(5, {{0, 3, 1}}),
      Perm(5, {{1, 3, 2}})
    }
  );

  PermGroup pg3(
    {
      Perm(5, {{2, 3, 0}})
    }
  );

  EXPECT_TRUE(pg1 == pg2 && !(pg1 != pg2))
    << "Can recognize permutation groups as equal.";

  EXPECT_TRUE(pg1 != pg3 && pg2 != pg3 && !(pg1 == pg3) && !(pg2 == pg3))
    << "Can recognize permutation groups as unequal.";
}

TEST(PermGroupTest, CanObtainDegree)
{
  PermGroup pg({Perm(10u)});
  EXPECT_EQ(10u, pg.degree())
    << "Permutation group degree set correctly.";
}

TEST(PermGroupTest, CanObtainOrder)
{
  PermGroup id(10u, {});
  EXPECT_EQ(1u, id.order())
    << "Order set correctly for trivial permutation group.";

  for (unsigned i = 1u; i <= 10u; ++i) {
    EXPECT_EQ(util::factorial(i), PermGroup::symmetric(i).order())
      << "Order set correctly for symmetric group S" << i;
  }

  for (unsigned i = 1u; i <= 10u; ++i) {
    EXPECT_EQ(i, PermGroup::cyclic(i).order())
      << "Order set correctly for cyclic group Z" << i;
  }

  for (unsigned i = 2u; i <= 10u; i += 2u) {
    EXPECT_EQ(i, PermGroup::dihedral(i).order())
      << "Order set correctly for dihedral group D" << i;
  }
}

TEST(PermGroupTest, CanCheckForSymmetricGroup)
{
  for (unsigned i = 1u; i < 10; ++i) {
    EXPECT_TRUE(PermGroup::symmetric(i).is_symmetric())
      << "Symmetric group correctly identified as such";
  }
}

TEST(PermGroupTest, CanDetermineTransitivity)
{
  PermGroup transitive_group(
    {
      Perm(9, {{0, 1}}),
      Perm(9, {{1, 2}}),
      Perm(9, {{2, 3, 4}}),
      Perm(9, {{4, 5, 6, 7, 8}})
    }
  );

  EXPECT_TRUE(transitive_group.is_transitive())
    << "Transitive group correctly identified as such.";

  PermGroup non_transitive_group(
    {
      Perm(14, {{0, 1}}),
      Perm(14, {{1, 2}}),
      Perm(14, {{3, 4}}),
      Perm(14, {{4, 5}}),
      Perm(14, {{6, 7}}),
      Perm(14, {{7, 8}}),
      Perm(14, {{11, 12}, {0, 3}, {1, 4}, {2, 5}}),
      Perm(14, {{12, 13}, {3, 6}, {4, 7}, {5, 8}})
    }
  );

  EXPECT_FALSE(non_transitive_group.is_transitive())
    << "Non-transitive group correctly identified as such.";
}

TEST(PermGroupTest, CanTestMembership)
{
  PermGroup a4(verified_perm_group(A4));

  std::vector<Perm> expected_members {
    Perm(4),
    Perm(4, {{1, 2, 3}}),
    Perm(4, {{1, 3, 2}}),
    Perm(4, {{0, 1}, {2, 3}}),
    Perm(4, {{0, 1, 2}}),
    Perm(4, {{0, 1, 3}}),
    Perm(4, {{0, 2, 1}}),
    Perm(4, {{0, 2, 3}}),
    Perm(4, {{0, 2}, {1, 3}}),
    Perm(4, {{0, 3, 1}}),
    Perm(4, {{0, 3, 2}}),
    Perm(4, {{0, 3}, {1, 2}})
  };

  std::vector<Perm> expected_non_members {
    Perm(4, {{2, 3}}),
    Perm(4, {{1, 2}}),
    Perm(4, {{1, 3}}),
    Perm(4, {{0, 1}}),
    Perm(4, {{0, 1, 2, 3}}),
    Perm(4, {{0, 1, 3, 2}}),
    Perm(4, {{0, 2, 3, 1}}),
    Perm(4, {{0, 2}}),
    Perm(4, {{0, 2, 1, 3}}),
    Perm(4, {{0, 3, 2, 1}}),
    Perm(4, {{0, 3}}),
    Perm(4, {{0, 3, 1, 2}})
  };

  for (Perm const &perm : expected_members) {
    EXPECT_TRUE(a4.contains_element(perm))
      << "Membership test correctly identifies group member " << perm;
  }

  for (Perm const &perm : expected_non_members) {
    EXPECT_FALSE(a4.contains_element(perm))
      << "Membership test correctly rejects non group member " << perm;
  }
}

TEST(PermGroupTest, CanGenerateRandomElement)
{
  PermGroup a4(verified_perm_group(A4));

  for (unsigned i = 0u; i < 1000u; ++i) {
    EXPECT_TRUE(a4.contains_element(a4.random_element()))
      << "Randomly generated group element is actually inside group.";
  }
}

TEST(PermGroupTest, CanIterateTrivialGroup)
{
  PermGroup id = PermGroup(4, {});

  std::vector<Perm> actual_members1;
  for (Perm const &p : id)
    actual_members1.push_back(p);

  ASSERT_EQ(1u, actual_members1.size())
    << "Iterating trivial permutation group yields one element (ranged for).";

  EXPECT_TRUE(perm_equal({0, 1, 2, 3}, actual_members1[0]))
    << "Iterating trivial permutation group yields identity permutation (ranged for)..";

  std::vector<Perm> actual_members2;
  for (PermGroup::const_iterator it = id.begin(); it != id.end(); ++it)
    actual_members2.push_back(*it);

  ASSERT_EQ(1u, actual_members2.size())
    << "Iterating trivial permutation group yields one element (explicit iterator).";

  EXPECT_TRUE(perm_equal({0, 1, 2, 3}, actual_members2[0]))
    << "Iterating trivial permutation group yields identity permutation (explicit iterator).";
}

TEST(PermGroupTest, CanIterateSimplestNonTrivialGroup)
{
  PermGroup pg = PermGroup(4, {Perm(4, {{0, 1}})});

  std::vector<Perm> expected_members {
    Perm(4),
    Perm(4, {{0, 1}})
  };

  std::vector<Perm> actual_members1;
  for (Perm const &p : pg)
    actual_members1.push_back(p);

  ASSERT_EQ(expected_members.size(), actual_members1.size())
    << "Iterating simplest non-trivial permutation group yields two element (ranged for).";

  EXPECT_THAT(actual_members1, UnorderedElementsAreArray(expected_members))
    << "Iterating simplest non-trivial permutation group yields correct permutation (ranged for)..";

  std::vector<Perm> actual_members2;
  for (PermGroup::const_iterator it = pg.begin(); it != pg.end(); ++it)
    actual_members2.push_back(*it);

  ASSERT_EQ(expected_members.size(), actual_members2.size())
    << "Iterating simplest non-trivial permutation group yields two element (explicit iterator).";

  EXPECT_THAT(actual_members2, UnorderedElementsAreArray(expected_members))
    << "Iterating simplest non-trivial permutation group yields correct permutation (explicit iterator).";
}

TEST(PermGroupTest, CanIterateElements)
{
  PermGroup a4(verified_perm_group(A4));

  std::vector<Perm> expected_members {
    Perm(4),
    Perm(4, {{1, 2, 3}}),
    Perm(4, {{1, 3, 2}}),
    Perm(4, {{0, 1}, {2, 3}}),
    Perm(4, {{0, 1, 2}}),
    Perm(4, {{0, 1, 3}}),
    Perm(4, {{0, 2, 1}}),
    Perm(4, {{0, 2, 3}}),
    Perm(4, {{0, 2}, {1, 3}}),
    Perm(4, {{0, 3, 1}}),
    Perm(4, {{0, 3, 2}}),
    Perm(4, {{0, 3}, {1, 2}})
  };

  std::vector<Perm> actual_members1;

  for (Perm const &p : a4)
    actual_members1.push_back(p);

  EXPECT_THAT(actual_members1, UnorderedElementsAreArray(expected_members))
    << "Iteration produces every element exactly once (ranged for).";

  std::vector<Perm> actual_members2;

  for (PermGroup::const_iterator it = a4.begin(); it != a4.end(); ++it) {
    EXPECT_EQ(4u, it->degree())
      << "Iterator dereferencing works correctly.";

    EXPECT_TRUE(it == it && it != a4.end())
      << "Iterator comparison works correctly.";

    actual_members2.push_back(*it);
  }

  EXPECT_THAT(actual_members2, UnorderedElementsAreArray(expected_members))
    << "Iteration produces every element exactly once (explicit iterator).";
}

class PermGroupConstructionMethodTest : public testing::TestWithParam<
  std::tuple<BSGSOptions::Construction, BSGSOptions::Transversals>> {};

TEST_P(PermGroupConstructionMethodTest, CanGenerateCorrectGroupElements)
{
  BSGSOptions bsgs_options;

  std::tie(bsgs_options.construction, bsgs_options.transversals) = GetParam();

  PermGroup groups[] = {
    PermGroup(BSGS(
      {
        Perm(4, {{1, 3}}),
        Perm(4, {{0, 1}, {2, 3}})
      },
      &bsgs_options)
    ),
    PermGroup(BSGS(
      {
        Perm(5, {{1, 3}, {2, 4}}),
        Perm(5, {{0, 1, 2, 4, 3}})
      },
      &bsgs_options)
    ),
    PermGroup(BSGS(
      {
        Perm(6, {{0, 1, 2, 3, 4, 5}})
      },
      &bsgs_options)
    ),
    PermGroup(BSGS(
      {
        Perm(7, {{1, 4}, {2, 5}, {3, 6}}),
        Perm(7, {{0, 1, 3, 2, 5, 6, 4}})
      },
      &bsgs_options)
    )
  };

  PermSet expected_elements[] = {
    {
      Perm(4, {{0, 1, 2, 3}}),
      Perm(4, {{0, 1}, {2, 3}}),
      Perm(4, {{0, 2}, {1, 3}}),
      Perm(4, {{0, 2}}),
      Perm(4, {{0, 3, 2, 1}}),
      Perm(4, {{0, 3}, {1, 2}}),
      Perm(4, {{1, 3}})
    },
    {
      Perm(5, {{1, 3}, {2, 4}}),
      Perm(5, {{0, 1}, {2, 3}}),
      Perm(5, {{0, 1, 2, 4, 3}}),
      Perm(5, {{0, 2}, {3, 4}}),
      Perm(5, {{0, 2, 3, 1, 4}}),
      Perm(5, {{0, 3, 4, 2, 1}}),
      Perm(5, {{0, 3}, {1, 4}}),
      Perm(5, {{0, 4}, {1, 2}}),
      Perm(5, {{0, 4, 1, 3, 2}})
    },
    {
      Perm(6, {{0, 1, 2, 3, 4, 5}}),
      Perm(6, {{0, 2, 4}, {1, 3, 5}}),
      Perm(6, {{0, 3}, {1, 4}, {2, 5}}),
      Perm(6, {{0, 4, 2}, {1, 5, 3}}),
      Perm(6, {{0, 5, 4, 3, 2, 1}})
    },
    {
      Perm(7, {{1, 4}, {2, 5}, {3, 6}}),
      Perm(7, {{0, 1}, {2, 6}, {3, 4}}),
      Perm(7, {{0, 1, 3, 2, 5, 6, 4}}),
      Perm(7, {{0, 2}, {1, 3}, {4, 5}}),
      Perm(7, {{0, 2, 4, 3, 6, 1, 5}}),
      Perm(7, {{0, 3}, {2, 4}, {5, 6}}),
      Perm(7, {{0, 3, 5, 4, 1, 2, 6}}),
      Perm(7, {{0, 4, 6, 5, 2, 3, 1}}),
      Perm(7, {{0, 4}, {1, 6}, {3, 5}}),
      Perm(7, {{0, 5}, {1, 2}, {4, 6}}),
      Perm(7, {{0, 5, 1, 6, 3, 4, 2}}),
      Perm(7, {{0, 6, 2, 1, 4, 5, 3}}),
      Perm(7, {{0, 6}, {1, 5}, {2, 3}})
    },
  };

  for (auto i = 0u; i < sizeof(groups) / sizeof(groups[0]); ++i) {
    EXPECT_TRUE(perm_group_equal(expected_elements[i], groups[i]))
      << "Group generated correctly";
  }
}

INSTANTIATE_TEST_SUITE_P(ConstructionMethods, PermGroupConstructionMethodTest,
  testing::Combine(
    testing::Values(BSGSOptions::Construction::SCHREIER_SIMS,
                    BSGSOptions::Construction::SCHREIER_SIMS_RANDOM),
    testing::Values(BSGSOptions::Transversals::EXPLICIT,
                    BSGSOptions::Transversals::SCHREIER_TREES)));
                    // TODO: SHALLOW_SCHREIER_TREES

TEST(PermGroupCombinationTest, CanConstructDirectProduct)
{
  std::vector<std::vector<PermGroup>> direct_products {
    {
      PermGroup(
        {
          Perm(3, {{0, 1}}),
          Perm(3, {{0, 1, 2}})
        }
      ),
      PermGroup(
        {
          Perm(3, {{0, 1}}),
          Perm(3, {{0, 1, 2}})
        }
      )
    }
  };

  PermSet expected_direct_products[] = {
    {
      Perm(6, {{0, 1, 2}, {3, 4, 5}}),
      Perm(6, {{0, 1, 2}, {3, 4}}),
      Perm(6, {{0, 1, 2}, {3, 5, 4}}),
      Perm(6, {{0, 1, 2}, {3, 5}}),
      Perm(6, {{0, 1, 2}, {4, 5}}),
      Perm(6, {{0, 1, 2}}),
      Perm(6, {{0, 1}, {3, 4, 5}}),
      Perm(6, {{0, 1}, {3, 4}}),
      Perm(6, {{0, 1}, {3, 5, 4}}),
      Perm(6, {{0, 1}, {3, 5}}),
      Perm(6, {{0, 1}, {4, 5}}),
      Perm(6, {{0, 1}}),
      Perm(6, {{0, 2, 1}, {3, 4, 5}}),
      Perm(6, {{0, 2, 1}, {3, 4}}),
      Perm(6, {{0, 2, 1}, {3, 5, 4}}),
      Perm(6, {{0, 2, 1}, {3, 5}}),
      Perm(6, {{0, 2, 1}, {4, 5}}),
      Perm(6, {{0, 2, 1}}),
      Perm(6, {{0, 2}, {3, 4, 5}}),
      Perm(6, {{0, 2}, {3, 4}}),
      Perm(6, {{0, 2}, {3, 5, 4}}),
      Perm(6, {{0, 2}, {3, 5}}),
      Perm(6, {{0, 2}, {4, 5}}),
      Perm(6, {{0, 2}}),
      Perm(6, {{1, 2}, {3, 4, 5}}),
      Perm(6, {{1, 2}, {3, 4}}),
      Perm(6, {{1, 2}, {3, 5, 4}}),
      Perm(6, {{1, 2}, {3, 5}}),
      Perm(6, {{1, 2}, {4, 5}}),
      Perm(6, {{1, 2}}),
      Perm(6, {{3, 4, 5}}),
      Perm(6, {{3, 4}}),
      Perm(6, {{3, 5, 4}}),
      Perm(6, {{3, 5}}),
      Perm(6, {{4, 5}})
    }
  };

  for (auto i = 0u; i < direct_products.size(); ++i) {
    auto groups(direct_products[i]);

    EXPECT_TRUE(perm_group_equal(
      expected_direct_products[i],
      PermGroup::direct_product(groups.begin(), groups.end())))
        << "Direct product construction correct.";
  }
}

TEST(PermGroupCombinationTest, CanConstructWreathProduct)
{
  std::vector<std::pair<PermGroup, PermGroup>> wreath_products {
    {
      PermGroup(
        {
          Perm(5, {{0, 2, 1}}),
          Perm(5, {{3, 4}})
        }
      ),
      PermGroup(
        {
          Perm(5, {{0, 2, 1}, {3, 4}})
        }
      )
    },
    {
      PermGroup(
        {
          Perm(9, {{0, 3, 6}}),
          Perm(9, {{0, 4, 8}})
        }
      ),
      PermGroup(
        {
          Perm(3, {{0, 1, 2}})
        }
      )
    },
    {
      PermGroup(
        {
          Perm(5, {{1, 3}}),
          Perm(5, {{2, 4}})
        }
      ),
      PermGroup(
        {
          Perm(4, {{0, 1, 3}}),
          Perm(4, {{2, 1}}),
        }
      )
    }
  };

  PermGroup expected_wreath_products[] = {
    PermGroup(
      {
        Perm(25, {{0, 2, 1}}),
        Perm(25, {{3, 4}}),
        Perm(25, {{5, 7, 6}}),
        Perm(25, {{8, 9}}),
        Perm(25, {{10, 12, 11}}),
        Perm(25, {{13, 14}}),
        Perm(25, {{15, 17, 16}}),
        Perm(25, {{18, 19}}),
        Perm(25, {{20, 22, 21}}),
        Perm(25, {{23, 24}}),
        Perm(25, {{0, 10, 5}, {1, 11, 6},
                  {2, 12, 7}, {3, 13, 8},
                  {4, 14, 9}, {15, 20},
                  {16, 21}, {17, 22},
                  {18, 23}, {19, 24}})
      }
    ),
    PermGroup(
      {
        Perm(27, {{0, 3, 6}}),
        Perm(27, {{0, 4, 8}}),
        Perm(27, {{0, 9, 18}, {1, 10, 19},
                  {2, 11, 20}, {3, 12, 21},
                  {4, 13, 22}, {5, 14, 23},
                  {6, 15, 24}, {7, 16, 25},
                  {8, 17, 26}}),
        Perm(27, {{9, 12, 15}}),
        Perm(27, {{9, 13, 17}}),
        Perm(27, {{18, 21, 24}}),
        Perm(27, {{18, 22, 26}})
      }
    ),
    PermGroup(
      {
        Perm(20, {{1, 3}}),
        Perm(20, {{2, 4}}),
        Perm(20, {{6, 8}}),
        Perm(20, {{7, 9}}),
        Perm(20, {{11, 13}}),
        Perm(20, {{12, 14}}),
        Perm(20, {{16, 18}}),
        Perm(20, {{17, 19}}),
        Perm(20, {{0, 5, 15}, {1, 6, 16},
                  {2, 7, 17}, {3, 8, 18},
                  {4, 9, 19}}),
        Perm(20, {{5,10}, {6, 11},
                  {7, 12}, {8, 13},
                  {9, 14}}),
      }
    )
  };

  for (auto i = 0u; i < wreath_products.size(); ++i) {
    auto lhs(wreath_products[i].first);
    auto rhs(wreath_products[i].second);

    EXPECT_EQ(expected_wreath_products[i],
              PermGroup::wreath_product(lhs, rhs))
      << "Wreath product construction correct.";
  }
}

class DisjointSubgroupProductTest :
  public testing::TestWithParam<std::pair<bool, bool>> {};

TEST_P(DisjointSubgroupProductTest, CanFindDisjointSubgroupProduct)
{
  auto flags = GetParam();
  bool complete = std::get<0>(flags);
  bool disjoint_opt = std::get<1>(flags);

  std::vector<PermGroup> permgroups {
    PermGroup(
      {
        Perm(14, {{0, 1}}),
        Perm(14, {{1, 2}}),
        Perm(14, {{3, 4}}),
        Perm(14, {{4, 5}}),
        Perm(14, {{6, 7}}),
        Perm(14, {{7, 8}}),
        Perm(14, {{11, 12}, {0, 3}, {1, 4}, {2, 5}}),
        Perm(14, {{12, 13}, {3, 6}, {4, 7}, {5, 8}}),
        Perm(14, {{9, 10}})
      }
    )
  };

  if (complete) {
    permgroups.push_back(
      PermGroup(
        {
          Perm(21, {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}, {9, 10, 11}, {13, 14},
                    {16, 17}, {19, 20}}),
          Perm(21, {{1, 2}, {4, 5}, {7, 8}, {10, 11}, {12, 13, 14},
                    {15, 16, 17}, {18, 19, 20}})
        }
      )
    );
  }

  std::vector<std::vector<PermGroup>> expected_disjoint_subgroups {
    {
      PermGroup(
        {
          Perm(14, {{0, 1}}),
          Perm(14, {{1, 2}}),
          Perm(14, {{3, 4}}),
          Perm(14, {{4, 5}}),
          Perm(14, {{6, 7}}),
          Perm(14, {{7, 8}}),
          Perm(14, {{11, 12}, {0, 3}, {1, 4}, {2, 5}}),
          Perm(14, {{12, 13}, {3, 6}, {4, 7}, {5, 8}})
        }
      ),
      PermGroup(
        {
          Perm(14, {{9, 10}})
        }
      )
    }
  };

  if (complete) {
    expected_disjoint_subgroups.push_back(
      {
        PermGroup(
          {
            Perm(21, {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}, {9, 10, 11}}),
            Perm(21, {{0, 1}, {3, 4}, {6, 7}, {9, 10}})
          }
        ),
        PermGroup(
          {
            Perm(21, {{12, 13, 14}, {15, 16, 17}, {18, 19, 20}}),
            Perm(21, {{12, 13}, {15, 16}, {18, 19}})
          }
        ),
      }
    );
  }

  for (auto i = 0u; i < permgroups.size(); ++i) {
    auto disjoint_subgroups =
      permgroups[i].disjoint_decomposition(complete, disjoint_opt);

    EXPECT_THAT(disjoint_subgroups,
                UnorderedElementsAreArray(expected_disjoint_subgroups[i]))
      << "Disjoint subgroup product decomposition generated correctly.";
  }
}

INSTANTIATE_TEST_SUITE_P(DisjointSubgroupProductVariants,
  DisjointSubgroupProductTest, testing::Values(std::make_pair(false, false),
                                               std::make_pair(true, false),
                                               std::make_pair(true, true)));

//TEST(DISABLED_WreathProductTest, CanFindWreathProduct)
//{
//  PermGroup pg(12,
//    {
//      Perm(12, {{1, 2}}),
//      Perm(12, {{2, 3}}),
//      Perm(12, {{4, 5}}),
//      Perm(12, {{5, 6}}),
//      Perm(12, {{7, 8}}),
//      Perm(12, {{8, 9}}),
//      Perm(12, {{1, 4}, {2, 5}, {3, 6}, {10, 11}}),
//      Perm(12, {{4, 7}, {5, 8}, {6, 9}, {11, 12}})
//    }
//  );
//
//  std::vector<PermGroup> decomp = pg.wreath_decomposition();
//
//  ASSERT_EQ(4u, decomp.size())
//    << "Wreath product decomposition found.";
//
//  PermGroup sigma_k(12,
//    {
//      Perm(12, {{1, 4}, {2, 5}, {3, 6}, {10, 11}}),
//      Perm(12, {{4, 7}, {5, 8}, {6, 9}, {11, 12}})
//    }
//  );
//
//  EXPECT_EQ(sigma_k, decomp[0])
//    << "Block permuter monomorphism image generated correctly.";
//
//  std::vector<PermGroup> sigma_hs {
//    PermGroup(12,
//      {
//        Perm(12, {{1, 2}}),
//        Perm(12, {{2, 3}})
//      }
//    ),
//    PermGroup(12,
//      {
//        Perm(12, {{4, 5}}),
//        Perm(12, {{5, 6}})
//      }
//    ),
//    PermGroup(12,
//      {
//        Perm(12, {{7, 8}}),
//        Perm(12, {{8, 9}})
//      }
//    )
//  };
//
//  std::vector<PermGroup> tmp(decomp.begin() + 1, decomp.end());
//  EXPECT_THAT(tmp, UnorderedElementsAreArray(sigma_hs))
//    << "Permutation representations of block actions generated correctly.";
//}
//
//TEST(SpecialPermGroupTest, CanConstructSymmetricGroup)
//{
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(S1),
//                               PermGroup::symmetric(1)))
//    << "Can construct symmetric group S_1.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(S2),
//                               PermGroup::symmetric(2)))
//    << "Can construct symmetric group S_2.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(S3),
//                               PermGroup::symmetric(3)))
//    << "Can construct symmetric group S_3.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(S4),
//                               PermGroup::symmetric(4)))
//    << "Can construct symmetric group S_4.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(S5),
//                               PermGroup::symmetric(5)))
//    << "Can construct symmetric group S_5.";
//}
//
//TEST(SpecialPermGroupTest, CanConstructCyclicGroup)
//{
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(C1),
//                               PermGroup::cyclic(1)))
//    << "Can construct cyclic group C_1.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(C2),
//                               PermGroup::cyclic(2)))
//    << "Can construct cyclic group C_2.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(C3),
//                               PermGroup::cyclic(3)))
//    << "Can construct cyclic group C_3.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(C4),
//                               PermGroup::cyclic(4)))
//    << "Can construct cyclic group C_4.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(C5),
//                               PermGroup::cyclic(5)))
//    << "Can construct cyclic group C_5.";
//}
//
//TEST(SpecialPermGroupTest, CanConstructDihedralGroup)
//{
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(D2),
//                               PermGroup::dihedral(2)))
//    << "Can construct dihedral group D_2.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(D4),
//                               PermGroup::dihedral(4)))
//    << "Can construct dihedral group D_4.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(D6),
//                               PermGroup::dihedral(6)))
//    << "Can construct dihedral group D_6.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(D8),
//                               PermGroup::dihedral(8)))
//    << "Can construct dihedral group D_8.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(D10),
//                               PermGroup::dihedral(10)))
//    << "Can construct dihedral group D_10.";
//
//  EXPECT_TRUE(perm_group_equal(verified_perm_group(D12),
//                               PermGroup::dihedral(12)))
//    << "Can construct dihedral group D_12.";
//}
