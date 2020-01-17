#include <algorithm>
#include <vector>

#include "gmock/gmock.h"

#include "perm.h"
#include "perm_set.h"
#include "pr_randomizer.h"

#include "test_main.cc"

#define RANDOMIZER_RUNS 10000
#define RANDOMIZER_EPS_REL 5

using cgtl::Perm;
using cgtl::PermSet;
using cgtl::PrRandomizer;

class PRRandomizerTest : public testing::Test
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
      EXPECT_THAT(pr_expected[i], testing::Contains(pr_randomizers[i].next()))
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

TEST_F(PRRandomizerTest, CanTestForAltSym)
{
  auto symmetric_generators = [](unsigned n) {
    if (n == 1u)
      return PermSet {Perm(1u)};

    std::vector<unsigned> chain;
    for (unsigned i = 1u; i <= n; ++i)
      chain.push_back(i);

    return PermSet {Perm(n, {{1, 2}}), Perm(n, {chain})};
  };

  auto alternating_generators = [](unsigned n) {
    PermSet gens;
    for (unsigned i = 3u; i <= n; ++i)
      gens.emplace(Perm(n, {{1, 2, i}}));

    return gens;
  };

  for (unsigned i = 8u; i < 20u; ++i) {
    auto gens_symmetric(symmetric_generators(i));
    PrRandomizer pr_symmetric(gens_symmetric);

    auto gens_alternating(alternating_generators(i));
    PrRandomizer pr_alternating(gens_alternating);

    EXPECT_TRUE(pr_symmetric.test_symmetric())
      << "Can identify symmetric generating sets.";

    EXPECT_FALSE(pr_alternating.test_symmetric())
      << "Can identify non-symmetric generating sets.";

    EXPECT_TRUE(pr_alternating.test_alternating())
      << "Can identify alternating generating sets.";

    EXPECT_FALSE(pr_symmetric.test_alternating())
      << "Can identify non-alternating generating sets.";
  }
}
