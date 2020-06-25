#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

#include "gmock/gmock.h"

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

TEST(PermGroupTest, CanComparePermGroups)
{
  PermGroup pg1(5,
    {
      Perm(5, {{1, 2}, {3, 4}}),
      Perm(5, {{1, 4, 2}})
    }
  );

  PermGroup pg2(5,
    {
      Perm(5, {{1, 2}, {3, 4}}),
      Perm(5, {{1, 4, 2}}),
      Perm(5, {{2, 4, 3}})
    }
  );

  PermGroup pg3(5,
    {
      Perm(5, {{3, 4, 1}})
    }
  );

  EXPECT_TRUE(pg1 == pg2 && !(pg1 != pg2))
    << "Can recognize permutation groups as equal.";

  EXPECT_TRUE(pg1 != pg3 && pg2 != pg3 && !(pg1 == pg3) && !(pg2 == pg3))
    << "Can recognize permutation groups as unequal.";
}

TEST(PermGroupTest, CanObtainDegree)
{
  PermGroup pg(10u, {Perm(10u)});
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

  for (unsigned i = 3u; i <= 10u; ++i) {
    EXPECT_EQ(util::factorial(i) / 2, PermGroup::alternating(i).order())
      << "Order set correctly for alternating group A" << i;
  }
}

TEST(PermGroupTest, CanCheckForSymmetricGroup)
{
  for (unsigned i = 1u; i < 10; ++i) {
    EXPECT_TRUE(PermGroup::symmetric(i).is_symmetric())
      << "Symmetric group correctly identified as such";
  }
}

TEST(PermGroupTest, CanCheckForAlternatingGroup)
{
  for (unsigned i = 3u; i < 10; ++i) {
    EXPECT_TRUE(PermGroup::alternating(i).is_alternating())
      << "Alternating group correctly identified as such";
  }
}

TEST(PermGroupTest, CanDetermineTrasitivity)
{
  PermGroup transitive_group(9,
    {
      Perm(9, {{1, 2}}),
      Perm(9, {{2, 3}}),
      Perm(9, {{3, 4, 5}}),
      Perm(9, {{5, 6, 7, 8, 9}})
    }
  );

  EXPECT_TRUE(transitive_group.is_transitive())
    << "Transitive group correctly identified as such.";

  PermGroup non_transitive_group(14,
    {
      Perm(14, {{1, 2}}),
      Perm(14, {{2, 3}}),
      Perm(14, {{4, 5}}),
      Perm(14, {{5, 6}}),
      Perm(14, {{7, 8}}),
      Perm(14, {{8, 9}}),
      Perm(14, {{12, 13}, {1, 4}, {2, 5}, {3, 6}}),
      Perm(14, {{13, 14}, {4, 7}, {5, 8}, {6, 9}})
    }
  );

  EXPECT_FALSE(non_transitive_group.is_transitive())
    << "Non-transitive group correctly identified as such.";
}

TEST(PermGroupTest, CanTestMembership)
{
  PermGroup a4 = PermGroup::alternating(4);

  std::vector<Perm> expected_members {
    Perm(4),
    Perm(4, {{2, 3, 4}}),
    Perm(4, {{2, 4, 3}}),
    Perm(4, {{1, 2}, {3, 4}}),
    Perm(4, {{1, 2, 3}}),
    Perm(4, {{1, 2, 4}}),
    Perm(4, {{1, 3, 2}}),
    Perm(4, {{1, 3, 4}}),
    Perm(4, {{1, 3}, {2, 4}}),
    Perm(4, {{1, 4, 2}}),
    Perm(4, {{1, 4, 3}}),
    Perm(4, {{1, 4}, {2, 3}})
  };

  std::vector<Perm> expected_non_members {
    Perm(4, {{3, 4}}),
    Perm(4, {{2, 3}}),
    Perm(4, {{2, 4}}),
    Perm(4, {{1, 2}}),
    Perm(4, {{1, 2, 3, 4}}),
    Perm(4, {{1, 2, 4, 3}}),
    Perm(4, {{1, 3, 4, 2}}),
    Perm(4, {{1, 3}}),
    Perm(4, {{1, 3, 2, 4}}),
    Perm(4, {{1, 4, 3, 2}}),
    Perm(4, {{1, 4}}),
    Perm(4, {{1, 4, 2, 3}})
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
  PermGroup a4 = PermGroup::alternating(4);

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

  EXPECT_TRUE(perm_equal({1, 2, 3, 4}, actual_members1[0]))
    << "Iterating trivial permutation group yields identity permutation (ranged for)..";

  std::vector<Perm> actual_members2;
  for (PermGroup::const_iterator it = id.begin(); it != id.end(); ++it)
    actual_members2.push_back(*it);

  ASSERT_EQ(1u, actual_members2.size())
    << "Iterating trivial permutation group yields one element (explicit iterator).";

  EXPECT_TRUE(perm_equal({1, 2, 3, 4}, actual_members2[0]))
    << "Iterating trivial permutation group yields identity permutation (explicit iterator).";
}

TEST(PermGroupTest, CanIterateSimplestNonTrivialGroup)
{
  PermGroup pg = PermGroup(4, {Perm(4, {{1, 2}})});

  std::vector<Perm> expected_members {
    Perm(4),
    Perm(4, {{1, 2}})
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
  PermGroup a4 = PermGroup::alternating(4);

  std::vector<Perm> expected_members {
    Perm(4),
    Perm(4, {{2, 3, 4}}),
    Perm(4, {{2, 4, 3}}),
    Perm(4, {{1, 2}, {3, 4}}),
    Perm(4, {{1, 2, 3}}),
    Perm(4, {{1, 2, 4}}),
    Perm(4, {{1, 3, 2}}),
    Perm(4, {{1, 3, 4}}),
    Perm(4, {{1, 3}, {2, 4}}),
    Perm(4, {{1, 4, 2}}),
    Perm(4, {{1, 4, 3}}),
    Perm(4, {{1, 4}, {2, 3}})
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
    PermGroup(BSGS(4,
      {
        Perm(4, {{2, 4}}),
        Perm(4, {{1, 2}, {3, 4}})
      },
      &bsgs_options)
    ),
    PermGroup(BSGS(5,
      {
        Perm(5, {{2, 4}, {3, 5}}),
        Perm(5, {{1, 2, 3, 5, 4}})
      },
      &bsgs_options)
    ),
    PermGroup(BSGS(6,
      {
        Perm(6, {{1, 2, 3, 4, 5, 6}})
      },
      &bsgs_options)
    ),
    PermGroup(BSGS(7,
      {
        Perm(7, {{2, 5}, {3, 6}, {4, 7}}),
        Perm(7, {{1, 2, 4, 3, 6, 7, 5}})
      },
      &bsgs_options)
    )
  };

  PermSet expected_elements[] = {
    {
      Perm(4, {{1, 2, 3, 4}}),
      Perm(4, {{1, 2}, {3, 4}}),
      Perm(4, {{1, 3}, {2, 4}}),
      Perm(4, {{1, 3}}),
      Perm(4, {{1, 4, 3, 2}}),
      Perm(4, {{1, 4}, {2, 3}}),
      Perm(4, {{2, 4}})
    },
    {
      Perm(5, {{2, 4}, {3, 5}}),
      Perm(5, {{1, 2}, {3, 4}}),
      Perm(5, {{1, 2, 3, 5, 4}}),
      Perm(5, {{1, 3}, {4, 5}}),
      Perm(5, {{1, 3, 4, 2, 5}}),
      Perm(5, {{1, 4, 5, 3, 2}}),
      Perm(5, {{1, 4}, {2, 5}}),
      Perm(5, {{1, 5}, {2, 3}}),
      Perm(5, {{1, 5, 2, 4, 3}})
    },
    {
      Perm(6, {{1, 2, 3, 4, 5, 6}}),
      Perm(6, {{1, 3, 5}, {2, 4, 6}}),
      Perm(6, {{1, 4}, {2, 5}, {3, 6}}),
      Perm(6, {{1, 5, 3}, {2, 6, 4}}),
      Perm(6, {{1, 6, 5, 4, 3, 2}})
    },
    {
      Perm(7, {{2, 5}, {3, 6}, {4, 7}}),
      Perm(7, {{1, 2}, {3, 7}, {4, 5}}),
      Perm(7, {{1, 2, 4, 3, 6, 7, 5}}),
      Perm(7, {{1, 3}, {2, 4}, {5, 6}}),
      Perm(7, {{1, 3, 5, 4, 7, 2, 6}}),
      Perm(7, {{1, 4}, {3, 5}, {6, 7}}),
      Perm(7, {{1, 4, 6, 5, 2, 3, 7}}),
      Perm(7, {{1, 5, 7, 6, 3, 4, 2}}),
      Perm(7, {{1, 5}, {2, 7}, {4, 6}}),
      Perm(7, {{1, 6}, {2, 3}, {5, 7}}),
      Perm(7, {{1, 6, 2, 7, 4, 5, 3}}),
      Perm(7, {{1, 7, 3, 2, 5, 6, 4}}),
      Perm(7, {{1, 7}, {2, 6}, {3, 4}})
    },
  };

  for (auto i = 0u; i < sizeof(groups) / sizeof(groups[0]); ++i) {
    EXPECT_TRUE(perm_group_equal(expected_elements[i], groups[i]))
      << "Group generated correctly";
  }
}

INSTANTIATE_TEST_CASE_P(ConstructionMethods, PermGroupConstructionMethodTest,
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
      PermGroup(3,
        {
          Perm(3, {{1, 2}}),
          Perm(3, {{1, 2, 3}})
        }
      ),
      PermGroup(3,
        {
          Perm(3, {{1, 2}}),
          Perm(3, {{1, 2, 3}})
        }
      )
    }
  };

  PermSet expected_direct_products[] = {
    {
      Perm(6, {{1, 2, 3}, {4, 5, 6}}),
      Perm(6, {{1, 2, 3}, {4, 5}}),
      Perm(6, {{1, 2, 3}, {4, 6, 5}}),
      Perm(6, {{1, 2, 3}, {4, 6}}),
      Perm(6, {{1, 2, 3}, {5, 6}}),
      Perm(6, {{1, 2, 3}}),
      Perm(6, {{1, 2}, {4, 5, 6}}),
      Perm(6, {{1, 2}, {4, 5}}),
      Perm(6, {{1, 2}, {4, 6, 5}}),
      Perm(6, {{1, 2}, {4, 6}}),
      Perm(6, {{1, 2}, {5, 6}}),
      Perm(6, {{1, 2}}),
      Perm(6, {{1, 3, 2}, {4, 5, 6}}),
      Perm(6, {{1, 3, 2}, {4, 5}}),
      Perm(6, {{1, 3, 2}, {4, 6, 5}}),
      Perm(6, {{1, 3, 2}, {4, 6}}),
      Perm(6, {{1, 3, 2}, {5, 6}}),
      Perm(6, {{1, 3, 2}}),
      Perm(6, {{1, 3}, {4, 5, 6}}),
      Perm(6, {{1, 3}, {4, 5}}),
      Perm(6, {{1, 3}, {4, 6, 5}}),
      Perm(6, {{1, 3}, {4, 6}}),
      Perm(6, {{1, 3}, {5, 6}}),
      Perm(6, {{1, 3}}),
      Perm(6, {{2, 3}, {4, 5, 6}}),
      Perm(6, {{2, 3}, {4, 5}}),
      Perm(6, {{2, 3}, {4, 6, 5}}),
      Perm(6, {{2, 3}, {4, 6}}),
      Perm(6, {{2, 3}, {5, 6}}),
      Perm(6, {{2, 3}}),
      Perm(6, {{4, 5, 6}}),
      Perm(6, {{4, 5}}),
      Perm(6, {{4, 6, 5}}),
      Perm(6, {{4, 6}}),
      Perm(6, {{5, 6}})
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
      PermGroup(5,
        {
          Perm(5, {{1, 3, 2}}),
          Perm(5, {{4, 5}})
        }
      ),
      PermGroup(5,
        {
          Perm(5, {{1, 3, 2}, {4, 5}})
        }
      )
    },
    {
      PermGroup(9,
        {
          Perm(9, {{1, 4, 7}}),
          Perm(9, {{1, 5, 9}})
        }
      ),
      PermGroup(3,
        {
          Perm(3, {{1, 2, 3}})
        }
      )
    },
    {
      PermGroup(5,
        {
          Perm(5, {{2, 4}}),
          Perm(5, {{3, 5}})
        }
      ),
      PermGroup(4,
        {
          Perm(4, {{1, 2, 4}}),
          Perm(4, {{3, 2}}),
        }
      )
    }
  };

  PermGroup expected_wreath_products[] = {
    PermGroup(25,
      {
        Perm(25, {{1, 3, 2}}),
        Perm(25, {{4, 5}}),
        Perm(25, {{6, 8, 7}}),
        Perm(25, {{9, 10}}),
        Perm(25, {{11, 13, 12}}),
        Perm(25, {{14, 15}}),
        Perm(25, {{16, 18, 17}}),
        Perm(25, {{19, 20}}),
        Perm(25, {{21, 23, 22}}),
        Perm(25, {{24, 25}}),
        Perm(25, {{1, 11, 6}, {2, 12, 7},
                  {3, 13, 8}, {4, 14, 9},
                  {5, 15, 10}, {16, 21},
                  {17, 22}, {18, 23},
                  {19, 24}, {20, 25}})
      }
    ),
    PermGroup(15,
      {
        Perm(15, {{1, 2, 4}}),
        Perm(15, {{1, 3, 5}}),
        Perm(15, {{1, 6, 11}, {2, 7, 12},
                  {3, 8, 13}, {4, 9, 14},
                  {5, 10, 15}}),
        Perm(15, {{11, 12, 14}}),
        Perm(15, {{11, 13, 15}}),
        Perm(15, {{6, 7, 9}}),
        Perm(15, {{6, 8, 10}})
      }
    ),
    PermGroup(16,
      {
        Perm(16, {{1, 3}}),
        Perm(16, {{1, 5, 13}, {2, 6, 14},
                  {3, 7, 15}, {4, 8, 16}}),
        Perm(16, {{10, 12}}),
        Perm(16, {{13, 15}}),
        Perm(16, {{14, 16}}),
        Perm(16, {{2, 4}}),
        Perm(16, {{5, 7}}),
        Perm(16, {{5, 9}, {6, 10},
                  {7, 11}, {8, 12}}),
        Perm(16, {{6, 8}}),
        Perm(16, {{9, 11}})
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
    PermGroup(14,
      {
        Perm(14, {{1, 2}}),
        Perm(14, {{2, 3}}),
        Perm(14, {{4, 5}}),
        Perm(14, {{5, 6}}),
        Perm(14, {{7, 8}}),
        Perm(14, {{8, 9}}),
        Perm(14, {{12, 13}, {1, 4}, {2, 5}, {3, 6}}),
        Perm(14, {{13, 14}, {4, 7}, {5, 8}, {6, 9}}),
        Perm(14, {{10, 11}})
      }
    )
  };

  if (complete) {
    permgroups.push_back(
      PermGroup(21,
        {
          Perm(21, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}, {14, 15},
                    {17, 18}, {20, 21}}),
          Perm(21, {{2, 3}, {5, 6}, {8, 9}, {11, 12}, {13, 14, 15},
                    {16, 17, 18}, {19, 20, 21}})
        }
      )
    );
  }

  std::vector<std::vector<PermGroup>> expected_disjoint_subgroups {
    {
      PermGroup(14,
        {
          Perm(14, {{1, 2}}),
          Perm(14, {{2, 3}}),
          Perm(14, {{4, 5}}),
          Perm(14, {{5, 6}}),
          Perm(14, {{7, 8}}),
          Perm(14, {{8, 9}}),
          Perm(14, {{12, 13}, {1, 4}, {2, 5}, {3, 6}}),
          Perm(14, {{13, 14}, {4, 7}, {5, 8}, {6, 9}})
        }
      ),
      PermGroup(14,
        {
          Perm(14, {{10, 11}})
        }
      )
    }
  };

  if (complete) {
    expected_disjoint_subgroups.push_back(
      {
        PermGroup(21,
          {
            Perm(21, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}),
            Perm(21, {{1, 2}, {4, 5}, {7, 8}, {10, 11}})
          }
        ),
        PermGroup(21,
          {
            Perm(21, {{13, 14, 15}, {16, 17, 18}, {19, 20, 21}}),
            Perm(21, {{13, 14}, {16, 17}, {19, 20}})
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

INSTANTIATE_TEST_CASE_P(DisjointSubgroupProductVariants,
  DisjointSubgroupProductTest, testing::Values(std::make_pair(false, false),
                                               std::make_pair(true, false),
                                               std::make_pair(true, true)));

TEST(DISABLED_WreathProductTest, CanFindWreathProduct)
{
  PermGroup pg(12,
    {
      Perm(12, {{1, 2}}),
      Perm(12, {{2, 3}}),
      Perm(12, {{4, 5}}),
      Perm(12, {{5, 6}}),
      Perm(12, {{7, 8}}),
      Perm(12, {{8, 9}}),
      Perm(12, {{1, 4}, {2, 5}, {3, 6}, {10, 11}}),
      Perm(12, {{4, 7}, {5, 8}, {6, 9}, {11, 12}})
    }
  );

  std::vector<PermGroup> decomp = pg.wreath_decomposition();

  ASSERT_EQ(4u, decomp.size())
    << "Wreath product decomposition found.";

  PermGroup sigma_k(12,
    {
      Perm(12, {{1, 4}, {2, 5}, {3, 6}, {10, 11}}),
      Perm(12, {{4, 7}, {5, 8}, {6, 9}, {11, 12}})
    }
  );

  EXPECT_EQ(sigma_k, decomp[0])
    << "Block permuter monomorphism image generated correctly.";

  std::vector<PermGroup> sigma_hs {
    PermGroup(12,
      {
        Perm(12, {{1, 2}}),
        Perm(12, {{2, 3}})
      }
    ),
    PermGroup(12,
      {
        Perm(12, {{4, 5}}),
        Perm(12, {{5, 6}})
      }
    ),
    PermGroup(12,
      {
        Perm(12, {{7, 8}}),
        Perm(12, {{8, 9}})
      }
    )
  };

  std::vector<PermGroup> tmp(decomp.begin() + 1, decomp.end());
  EXPECT_THAT(tmp, UnorderedElementsAreArray(sigma_hs))
    << "Permutation representations of block actions generated correctly.";
}

TEST(SpecialPermGroupTest, CanConstructSymmetricGroup)
{
  EXPECT_TRUE(perm_group_equal(verified_perm_group(S1),
                               PermGroup::symmetric(1)))
    << "Can construct symmetric group S_1.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(S2),
                               PermGroup::symmetric(2)))
    << "Can construct symmetric group S_2.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(S3),
                               PermGroup::symmetric(3)))
    << "Can construct symmetric group S_3.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(S4),
                               PermGroup::symmetric(4)))
    << "Can construct symmetric group S_4.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(S5),
                               PermGroup::symmetric(5)))
    << "Can construct symmetric group S_5.";
}

TEST(SpecialPermGroupTest, CanConstructCyclicGroup)
{
  EXPECT_TRUE(perm_group_equal(verified_perm_group(C1),
                               PermGroup::cyclic(1)))
    << "Can construct cyclic group C_1.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(C2),
                               PermGroup::cyclic(2)))
    << "Can construct cyclic group C_2.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(C3),
                               PermGroup::cyclic(3)))
    << "Can construct cyclic group C_3.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(C4),
                               PermGroup::cyclic(4)))
    << "Can construct cyclic group C_4.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(C5),
                               PermGroup::cyclic(5)))
    << "Can construct cyclic group C_5.";
}

TEST(SpecialPermGroupTest, CanConstructAlternatingGroup)
{
  EXPECT_TRUE(perm_group_equal(verified_perm_group(A1),
                               PermGroup::alternating(1)))
    << "Can construct alternating group A_1.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(A2),
                               PermGroup::alternating(2)))
    << "Can construct alternating group A_2.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(A3),
                               PermGroup::alternating(3)))
    << "Can construct alternating group A_3.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(A4),
                               PermGroup::alternating(4)))
    << "Can construct alternating group A_4.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(A5),
                               PermGroup::alternating(5)))
    << "Can construct alternating group A_5.";
}

TEST(SpecialPermGroupTest, CanConstructDihedralGroup)
{
  EXPECT_TRUE(perm_group_equal(verified_perm_group(D2),
                               PermGroup::dihedral(2)))
    << "Can construct dihedral group D_2.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D4),
                               PermGroup::dihedral(4)))
    << "Can construct dihedral group D_4.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D6),
                               PermGroup::dihedral(6)))
    << "Can construct dihedral group D_6.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D8),
                               PermGroup::dihedral(8)))
    << "Can construct dihedral group D_8.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D10),
                               PermGroup::dihedral(10)))
    << "Can construct dihedral group D_10.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D12),
                               PermGroup::dihedral(12)))
    << "Can construct dihedral group D_12.";
}
