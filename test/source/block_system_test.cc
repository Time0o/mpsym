#include <sstream>
#include <string>
#include <vector>

#include "gmock/gmock.h"

#include "block_system.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"

#include "test_main.cc"

using cgtl::BlockSystem;
using cgtl::Perm;
using cgtl::PermGroup;
using cgtl::PermSet;

static std::string block_to_string(std::vector<unsigned> const &block)
{
   std::stringstream ss;
   ss << '{' << block[0];
   for (auto i = 1u; i < block.size(); ++i)
     ss << ", " << block[i];
   ss << '}';

   return ss.str();
}

static testing::AssertionResult block_system_equal(
  std::vector<std::vector<unsigned>> const &expected, BlockSystem const &bs)
{
  if (static_cast<unsigned>(expected.size()) != bs.size()) {
    auto res = testing::AssertionFailure();
    return res << "Expected block system of size " << expected.size()
               << " but got one of size " << bs.size();
  }

  std::vector<int> block_found(expected.size(), 0);

  for (auto const &block : bs) {
     bool is_expected_block = false;
     for (auto i = 0u; i < expected.size(); ++i) {
       if (block_found[i])
         continue;

       testing::Matcher<std::vector<unsigned> const &> matcher(
         testing::UnorderedElementsAreArray(expected[i]));

       testing::StringMatchResultListener dummy;
       if (matcher.MatchAndExplain(block, &dummy)) {
         is_expected_block = true;
         block_found[i] = 1;
       }
     }

     if (!is_expected_block) {
       auto res = testing::AssertionFailure();
       return res << "Block " << block_to_string(block)
                  << " matches no expected block";
     }
  }

  for (auto i = 0u; i < expected.size(); ++i) {
    if (!block_found[i]) {
       auto res = testing::AssertionFailure();
       return res << "No match for block " << block_to_string(expected[i])
                  << " (more might be unmatched)";
    }
  }

  return testing::AssertionSuccess();
}

TEST(BlockSystemTest, CanFindMinimalBlockSystem)
{
  PermSet generators {
    Perm(6, {{1, 2, 3, 4, 5, 6}}), Perm(6, {{2, 6}, {3, 5}})
  };

  BlockSystem bs = BlockSystem::minimal(generators, {1, 3});

  EXPECT_TRUE(block_system_equal({{1, 3, 5}, {2, 4, 6}},
                                 BlockSystem::minimal(generators, {1, 3})))
    << "Minimal blocksystem correctly determined.";
}

TEST(BlockSystemTest, CanFindAllNonTrivialBlockSystemsForTransientGroup)
{
  PermGroup pg(9,
    {
      Perm(9, {{1, 2}}),
      Perm(9, {{1, 3}}),
      Perm(9, {{1, 4}, {2, 5}, {3, 6}}),
      Perm(9, {{1, 7}, {2, 8}, {3, 9}}),
      Perm(9, {{2, 3}}),
      Perm(9, {{4, 5}}),
      Perm(9, {{4, 7}, {5, 8}, {6, 9}}),
      Perm(9, {{5, 6}}),
      Perm(9, {{7, 8}}),
      Perm(9, {{7, 9}}),
      Perm(9, {{8, 9}})
    }
  );

  ASSERT_TRUE(pg.is_transitive())
    << "Permutation group is actually transitive.";

  auto block_systems(BlockSystem::non_trivial(pg, true));

  ASSERT_EQ(1u, block_systems.size())
    << "Correct number of block systems found.";

  EXPECT_TRUE(block_system_equal({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
              block_systems[0]))
    << "Correct block systems determined.";
}

TEST(BlockSystemTest, CanFindAllNonTrivialBlockSystemsForNonTransientGroup)
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

  ASSERT_FALSE(pg.is_transitive())
    << "Permutation group is actually non-transitive.";

  auto block_systems(BlockSystem::non_trivial(pg));

  ASSERT_EQ(1u, block_systems.size())
    << "Correct number of block systems found.";

  EXPECT_TRUE(block_system_equal({{1, 2, 3, 10}, {4, 5, 6, 11}, {7, 8, 9, 12}},
              block_systems[0]))
    << "Correct block systems determined.";
}
