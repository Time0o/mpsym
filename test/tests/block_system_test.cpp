#include <sstream>
#include <string>
#include <vector>

#include "gmock/gmock.h"

#include "block_system.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"

#include "test_main.cpp"

using namespace mpsym;
using namespace mpsym::internal;

static std::string block_to_string(BlockSystem::Block const &block)
{
   std::stringstream ss;
   ss << '{' << block[0];
   for (auto i = 1u; i < block.size(); ++i)
     ss << ", " << block[i];
   ss << '}';

   return ss.str();
}

static testing::AssertionResult block_system_equal(
  std::vector<BlockSystem::Block> const &expected, BlockSystem const &bs)
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

       testing::Matcher<BlockSystem::Block const &> matcher(
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
  std::vector<PermSet> generators {
    {
       Perm(6, {{0, 1, 2, 3, 4, 5}}),
       Perm(6, {{1, 5}, {2, 4}})
    },
    {
      Perm(9, {{0, 2}}),
      Perm(9, {{0, 3}, {1, 4}, {2, 5}}),
      Perm(9, {{3, 5}}),
      Perm(9, {{3, 6}, {4, 7}, {5, 8}}),
      Perm(9, {{6, 7}}),
      Perm(9, {{7, 8}})
    }
  };

  std::vector<std::vector<unsigned>> initial_classes {
    {0, 2},
    {0, 7}
  };

  std::vector<std::vector<BlockSystem::Block>> expected_block_systems {
    {{0, 2, 4}, {1, 3, 5}},
    {{0, 1, 2, 3, 4, 5, 6, 7, 8}}
  };

  for (auto i = 0u; i < generators.size(); ++i) {
    BlockSystem bs(BlockSystem::minimal(generators[i], initial_classes[i]));

    EXPECT_TRUE(block_system_equal(expected_block_systems[i], bs))
      << "Minimal blocksystem correctly determined.";
  }
}

TEST(BlockSystemTest, CanFindAllNonTrivialBlockSystemsForTransientGroup)
{
  PermGroup pg(
    {
      Perm(9, {{0, 1}}),
      Perm(9, {{0, 2}}),
      Perm(9, {{0, 3}, {1, 4}, {2, 5}}),
      Perm(9, {{0, 6}, {1, 7}, {2, 8}}),
      Perm(9, {{1, 2}}),
      Perm(9, {{3, 4}}),
      Perm(9, {{3, 6}, {4, 7}, {5, 8}}),
      Perm(9, {{4, 5}}),
      Perm(9, {{6, 7}}),
      Perm(9, {{6, 8}}),
      Perm(9, {{7, 8}})
    }
  );

  ASSERT_TRUE(pg.is_transitive())
    << "Permutation group is actually transitive.";

  auto block_systems(BlockSystem::non_trivial(pg, true));

  ASSERT_EQ(1u, block_systems.size())
    << "Correct number of block systems found.";

  EXPECT_TRUE(block_system_equal({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}},
              block_systems[0]))
    << "Correct block systems determined.";
}

TEST(BlockSystemTest, CanFindAllNonTrivialBlockSystemsForNonTransientGroup)
{
  PermGroup pg(
    {
      Perm(12, {{0, 1}}),
      Perm(12, {{1, 2}}),
      Perm(12, {{3, 4}}),
      Perm(12, {{4, 5}}),
      Perm(12, {{6, 7}}),
      Perm(12, {{7, 8}}),
      Perm(12, {{0, 3}, {1, 4}, {2, 5}, {9, 10}}),
      Perm(12, {{3, 6}, {4, 7}, {5, 8}, {10, 11}})
    }
  );

  ASSERT_FALSE(pg.is_transitive())
    << "Permutation group is actually non-transitive.";

  auto block_systems(BlockSystem::non_trivial(pg));

  ASSERT_EQ(1u, block_systems.size())
    << "Correct number of block systems found.";

  EXPECT_TRUE(block_system_equal({{0, 1, 2, 9}, {3, 4, 5, 10}, {6, 7, 8, 11}},
              block_systems[0]))
    << "Correct block systems determined.";
}
