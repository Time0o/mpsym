#include <sstream>
#include <utility>
#include <vector>

#include "gmock/gmock.h"

#include "perm.h"
#include "perm_group.h"
#include "schreier_sims.h"
#include "test_utility.h"

#include "test_main.cc"

using cgtl::Perm;
using cgtl::PermGroup;

using testing::ElementsAre;
using testing::UnorderedElementsAre;
using testing::UnorderedElementsAreArray;

static unsigned factorial(unsigned x)
{
  unsigned result = 1u;
  while (x > 1u)
    result *= x--;

  return result;
}

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
    EXPECT_EQ(factorial(i), PermGroup::symmetric(i).order())
      << "Order set correctly for symmetric group S" << i;
  }

  for (unsigned i = 1u; i <= 10u; ++i) {
    EXPECT_EQ(i, PermGroup::cyclic(i).order())
      << "Order set correctly for cyclic group Z" << i;
  }

  for (unsigned i = 3u; i <= 10u; ++i) {
    EXPECT_EQ(factorial(i) / 2, PermGroup::alternating(i).order())
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
  for (PermGroup::const_iterator it = id.begin(); it != id.end(); it++)
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
  for (PermGroup::const_iterator it = pg.begin(); it != pg.end(); it++)
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

  for (PermGroup::const_iterator it = a4.begin(); it != a4.end(); it++) {
    EXPECT_EQ(4u, it->degree())
      << "Iterator dereferencing works correctly.";

    EXPECT_TRUE(it == it && it != a4.end())
      << "Iterator comparison works correctly.";

    actual_members2.push_back(*it);
  }

  EXPECT_THAT(actual_members2, UnorderedElementsAreArray(expected_members))
    << "Iteration produces every element exactly once (explicit iterator).";
}

class PermGroupConstructionMethodTest
  : public testing::TestWithParam<PermGroup::ConstructionMethod> {};

// TODO: test more groups
TEST_P(PermGroupConstructionMethodTest, CanGenerateCorrectGroupElements)
{
  PermGroup::ConstructionMethod method = GetParam();

  typedef std::vector<std::vector<std::vector<unsigned>>> elemset;

  std::vector<std::pair<unsigned, elemset>> groups {
    {4, {{{2, 4}}, {{1, 2}, {3, 4}}}}
  };

  std::vector<elemset> expected_elements {
    {
      {{1, 2, 3, 4}},
      {{1, 2}, {3, 4}},
      {{1, 3}, {2, 4}},
      {{1, 3}},
      {{1, 4, 3, 2}},
      {{1, 4}, {2, 3}},
      {{2, 4}}
    }
  };

  for (auto i = 0u; i < groups.size(); ++i) {
    unsigned degree = std::get<0>(groups[i]);

    std::vector<Perm> generators;
    for (auto const &perm : std::get<1>(groups[i]))
      generators.push_back(Perm(degree, perm));

    std::stringstream ss;
    ss << "Group generated correctly, generators are: ";
    for (auto j = 0u; j < generators.size(); ++j)
      ss << generators[j] << (j == generators.size() - 1u ? "" : ", ");

    EXPECT_TRUE(perm_group_equal(
      expected_elements[i], PermGroup(degree, generators, method)))
      << ss.str();
  }
}

INSTANTIATE_TEST_CASE_P(ConstructionMethods, PermGroupConstructionMethodTest,
  testing::Values(PermGroup::ConstructionMethod::SCHREIER_SIMS,
                  PermGroup::ConstructionMethod::SCHREIER_SIMS_RANDOM));
                  // TODO: AUTO

TEST(PermGroupCombinationTest, CanConstructDirectProduct)
{
  auto S3(PermGroup::symmetric(3));

  std::vector<Perm> S3_shifted_generators;
  for (auto const &gen : S3.bsgs().strong_generators)
    S3_shifted_generators.push_back(gen.shifted(3));

  PermGroup S3_shifted(6, S3_shifted_generators);

  auto S3xS3(PermGroup::direct_product(S3, S3_shifted));
  auto S3xS3_autoshift(PermGroup::direct_product(S3, S3, true));

  std::vector<std::pair<PermGroup, PermGroup>> direct_products {
    {S3xS3, S3xS3_autoshift}
  };

  std::vector<std::vector<std::vector<unsigned>>> expected_direct_products[] =  {
    {
      {{1, 2, 3}, {4, 5, 6}},
      {{1, 2, 3}, {4, 5}},
      {{1, 2, 3}, {4, 6, 5}},
      {{1, 2, 3}, {4, 6}},
      {{1, 2, 3}, {5, 6}},
      {{1, 2, 3}},
      {{1, 2}, {4, 5, 6}},
      {{1, 2}, {4, 5}},
      {{1, 2}, {4, 6, 5}},
      {{1, 2}, {4, 6}},
      {{1, 2}, {5, 6}},
      {{1, 2}},
      {{1, 3, 2}, {4, 5, 6}},
      {{1, 3, 2}, {4, 5}},
      {{1, 3, 2}, {4, 6, 5}},
      {{1, 3, 2}, {4, 6}},
      {{1, 3, 2}, {5, 6}},
      {{1, 3, 2}},
      {{1, 3}, {4, 5, 6}},
      {{1, 3}, {4, 5}},
      {{1, 3}, {4, 6, 5}},
      {{1, 3}, {4, 6}},
      {{1, 3}, {5, 6}},
      {{1, 3}},
      {{2, 3}, {4, 5, 6}},
      {{2, 3}, {4, 5}},
      {{2, 3}, {4, 6, 5}},
      {{2, 3}, {4, 6}},
      {{2, 3}, {5, 6}},
      {{2, 3}},
      {{4, 5, 6}},
      {{4, 5}},
      {{4, 6, 5}},
      {{4, 6}},
      {{5, 6}},
    }
  };

  for (auto i = 0u; i < direct_products.size(); ++i) {
    EXPECT_TRUE(perm_group_equal(
      expected_direct_products[i], direct_products[i].first))
        << "Direct product construction correct (no autoshift).";

    EXPECT_TRUE(perm_group_equal(
      expected_direct_products[i], direct_products[i].second))
        << "Direct product construction correct (autoshift).";
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

TEST(WreathProductTest, CanFindWreathProduct)
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
  EXPECT_TRUE(perm_group_equal(verified_perm_group(D1),
                               PermGroup::dihedral(1)))
    << "Can construct dihedral group D_1.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D2),
                               PermGroup::dihedral(2)))
    << "Can construct dihedral group D_2.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D3),
                               PermGroup::dihedral(3)))
    << "Can construct dihedral group D_3.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D4),
                               PermGroup::dihedral(4)))
    << "Can construct dihedral group D_4.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D5),
                               PermGroup::dihedral(5)))
    << "Can construct dihedral group D_5.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D6),
                               PermGroup::dihedral(6)))
    << "Can construct dihedral group D_6.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D7),
                               PermGroup::dihedral(7)))
    << "Can construct dihedral group D_7.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D8),
                               PermGroup::dihedral(8)))
    << "Can construct dihedral group D_8.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D9),
                               PermGroup::dihedral(9)))
    << "Can construct dihedral group D_9.";

  EXPECT_TRUE(perm_group_equal(verified_perm_group(D10),
                               PermGroup::dihedral(10)))
    << "Can construct dihedral group D_10.";
}

TEST(SpecialPermGroupTest, CanConstructSymmetricGroupWithSupport)
{
  std::vector<PermGroup> symmetric_groups {
    PermGroup::symmetric({6, 9}),
    PermGroup::symmetric({7, 2, 4})
  };

  std::vector<std::vector<std::vector<unsigned>>> expected_elements[] = {
    {
      {{6, 9}}
    },
    {
      {{7, 2, 4}},
      {{7, 2}},
      {{7, 4, 2}},
      {{7, 4}},
      {{2, 4}}
    }
  };

  for (auto i = 0u; i < symmetric_groups.size(); ++i) {
    EXPECT_TRUE(perm_group_equal(expected_elements[i], symmetric_groups[i]))
      << "Symmetric group constructed for specific support has correct elements.";
  }
}

TEST(SpecialPermGroupTest, CanConstructCyclicGroupWithSupport)
{
  std::vector<PermGroup> cyclic_groups {
    PermGroup::cyclic({6, 9}),
    PermGroup::cyclic({7, 2, 4}),
    PermGroup::cyclic({1, 8, 4, 5})
  };

  std::vector<std::vector<std::vector<unsigned>>> expected_elements[] = {
    {
      {{6, 9}}
    },
    {
      {{7, 2, 4}},
      {{7, 4, 2}}
    },
    {
      {{1, 8, 4, 5}},
      {{1, 4}, {8, 5}},
      {{1, 5, 4, 8}}
    }
  };

  for (auto i = 0u; i < cyclic_groups.size(); ++i) {
    EXPECT_TRUE(perm_group_equal(expected_elements[i], cyclic_groups[i]))
      << "Cyclic group constructed for specific support has correct elements.";
  }
}

TEST(SpecialPermGroupTest, CanConstructAlternatingGroupWithSupport)
{
  std::vector<PermGroup> alternating_groups {
    PermGroup::alternating({7, 2, 4}),
    PermGroup::alternating({1, 8, 4, 5})
  };

  std::vector<std::vector<std::vector<unsigned>>> expected_elements[] = {
    {
      {{7, 2, 4}},
      {{7, 4, 2}}
    },
    {
      {{1, 8, 4}},
      {{1, 8, 5}},
      {{1, 8}, {4, 5}},
      {{1, 4, 8}},
      {{1, 4, 5}},
      {{1, 4}, {8, 5}},
      {{1, 5, 8}},
      {{1, 5, 4}},
      {{1, 5}, {8, 4}},
      {{8, 4, 5}},
      {{8, 5, 4}}
    }
  };

  for (auto i = 0u; i < alternating_groups.size(); ++i) {
    EXPECT_TRUE(perm_group_equal(expected_elements[i], alternating_groups[i]))
      << "Alternating group constructed for specific support has correct elements.";
  }
}

TEST(SpecialPermGroupTest, CanConstructDihedralGroupWithSupport)
{
  std::vector<PermGroup> dihedral_groups {
    PermGroup::dihedral({7, 2, 4}),
    PermGroup::dihedral({1, 8, 4, 5})
  };

  std::vector<std::vector<std::vector<unsigned>>> expected_elements[] = {
    {
      {{7, 2, 4}},
      {{7, 2}},
      {{7, 4, 2}},
      {{7, 4}},
      {{2, 4}}
    },
    {
      {{1, 8, 4, 5}},
      {{1, 8}, {4, 5}},
      {{1, 4}, {8, 5}},
      {{1, 4}},
      {{1, 5, 4, 8}},
      {{1, 5}, {8, 4}},
      {{8, 5}}
    }
  };

  for (auto i = 0u; i < dihedral_groups.size(); ++i) {
    EXPECT_TRUE(perm_group_equal(expected_elements[i], dihedral_groups[i]))
      << "Dihedral group constructed for specific support has correct elements.";
  }
}
