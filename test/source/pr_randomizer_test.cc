#include <algorithm>
#include <vector>

#include "gmock/gmock.h"
#include "pr_randomizer.h"

#include "test_main.cc"

#define RANDOMIZER_RUNS 10000
#define RANDOMIZER_EPS_REL 5

using cgtl::Perm;
using cgtl::PrRandomizer;

class PRRandomizerTest : public ::testing::Test
{
protected:
  std::vector<PrRandomizer> pr_randomizers {
    PrRandomizer({Perm(4, {{2, 4}}), Perm(4, {{1, 2}, {3, 4}})})
  };

  std::vector<std::vector<Perm>> pr_expected {
    std::vector<Perm> {
      Perm(4),
      Perm(4, {{1, 2, 3, 4}}),
      Perm(4, {{1, 3}, {2, 4}}),
      Perm(4, {{1, 4, 3, 2}}),
      Perm(4, {{1, 4}, {2, 3}}),
      Perm(4, {{1, 2}, {3, 4}}),
      Perm(4, {{1, 3}}),
      Perm(4, {{2, 4}})
    }
  };
};

TEST_F(PRRandomizerTest, CanContructRandomGroupMembers)
{
  for (auto i = 0u; i < pr_randomizers.size(); ++i) {
    for (int j = 0; j < RANDOMIZER_RUNS; ++j) {
      EXPECT_THAT(pr_expected[i], ::testing::Contains(pr_randomizers[i].next()))
        << "Product replacement randomizer produces group members.";
    }
  }
}

TEST_F(PRRandomizerTest, DistributionApproximatelyUniform)
{
  for (auto i = 0u; i < pr_randomizers.size(); ++i) {
    std::vector<unsigned> counts(pr_expected[i].size(), 0);

    for (int j = 0; j < RANDOMIZER_RUNS; ++j) {
      auto pos = std::find(pr_expected[i].begin(), pr_expected[i].end(),
                           pr_randomizers[i].next());

      ASSERT_NE(pr_expected[i].end(), pos)
        << "Generated element in group.";

      ++counts[pos - pr_expected[i].begin()];
    }

    unsigned expected_mean = RANDOMIZER_RUNS / pr_expected[i].size();
    unsigned allowed_delta = expected_mean / RANDOMIZER_EPS_REL;

    for (auto k = 0u; k < pr_expected[i].size(); ++k) {
      unsigned c = counts[k];

      unsigned min = expected_mean - allowed_delta;
      unsigned max = expected_mean + allowed_delta;

      EXPECT_TRUE(c >= min && c <= max)
        << "Value distribution approximately uniform (element "
        << pr_expected[i][k] << " occurred " << c << "/" << RANDOMIZER_RUNS
        << " times but should be in range [" << min << ", " << max << "])";
    }
  }
}
