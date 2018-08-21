#include <sstream>
#include <string>
#include <vector>

#include "block_system.h"
#include "gmock/gmock.h"
#include "perm.h"

#include "test_main.cc"

using cgtl::BlockSystem;
using cgtl::Perm;

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
  std::vector<Perm> generators {
    Perm(6, {{1, 2, 3, 4, 5, 6}}), Perm(6, {{2, 6}, {3, 5}})
  };

  BlockSystem bs = BlockSystem::minimal(generators, {1, 3});

  EXPECT_TRUE(block_system_equal({{1, 3, 5}, {2, 4, 6}},
                                 BlockSystem::minimal(generators, {1, 3})))
    << "Minimal blocksystem correctly determined.";
}
