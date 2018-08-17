#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>

#include "arch_graph.h"
#include "gmock/gmock.h"
#include "perm.h"
#include "perm_group.h"
#include "test_utility.h"

#include "test_main.cc"

using cgtl::ArchGraphSystem;
using cgtl::ArchGraph;
using cgtl::ArchGraphCluster;
using cgtl::Perm;
using cgtl::PermGroup;
using cgtl::TaskMapping;

using testing::UnorderedElementsAreArray;

typedef std::vector<std::vector<unsigned>> orbit;

static void expect_mapping_generates_orbits(
  ArchGraphSystem const *ag, std::vector<orbit> expected_orbits,
  ArchGraphSystem::MappingVariant mapping_variant)
{
  std::vector<TaskMapping> task_mappings;
  for (auto j = 0u; j < ag->num_processors(); ++j) {
    for (auto k = 0u; k < ag->num_processors(); ++k)
      task_mappings.push_back(ag->mapping({j, k}, 0u, mapping_variant));
  }

  for (auto const &tm1 : task_mappings) {
    std::vector<unsigned> mapping1(tm1.mapping());

    std::stringstream ss; ss << "{ ";
    for (auto j = 0u; j < mapping1.size(); ++j)
      ss << mapping1[j] << (j == mapping1.size() - 1u ? " }" : ", ");

    std::vector<std::vector<unsigned>> equivalent_assignments;
    for (auto const &tm2 : task_mappings) {
      if (tm1.equivalence_class() == tm2.equivalence_class())
        equivalent_assignments.push_back(tm2.mapping());
    }

    bool found_orbit = false;
    for (auto const &orbit : expected_orbits) {
      if (std::find(orbit.begin(), orbit.end(), tm1.mapping()) != orbit.end()) {

        EXPECT_THAT(equivalent_assignments, UnorderedElementsAreArray(orbit))
          << "Equivalent task mappings for " << ss.str() << " correct.";

        found_orbit = true;
        break;
      }
    }

    EXPECT_TRUE(found_orbit) << "Task mapping " << ss.str()
                             << " present in some orbit.";
  }
}

template<typename T>
class ArchGraphTestBase : public T
{
protected:
  void SetUp() {
    ag_nocol = construct_nocol();
    ag_vcol = construct_vcol();
    ag_ecol = construct_ecol();
    ag_tcol = construct_tcol();

    ag_nocol.todot(resource_path("ag_nocol.dot"));
    ag_vcol.todot(resource_path("ag_vcol.dot"));
    ag_ecol.todot(resource_path("ag_ecol.dot"));
    ag_tcol.todot(resource_path("ag_tcol.dot"));
  }

  ArchGraph ag_nocol;
  ArchGraph ag_vcol;
  ArchGraph ag_ecol;
  ArchGraph ag_tcol;

private:
  ArchGraph construct_nocol() {
    /*
     * 1 -- 1 -- 2  P -- C -- P
     * |         |  |         |
     * 4         2  C         C
     * |         |  |         |
     * 4 -- 3 -- 3  P -- C -- P
     */
    ArchGraph ag;

    auto p = ag.new_processor_type("P");
    auto c = ag.new_channel_type("C");

    auto pe1 = ag.add_processor(p);
    auto pe2 = ag.add_processor(p);
    auto pe3 = ag.add_processor(p);
    auto pe4 = ag.add_processor(p);

    ag.add_channel(pe1, pe2, c);
    ag.add_channel(pe2, pe3, c);
    ag.add_channel(pe3, pe4, c);
    ag.add_channel(pe4, pe1, c);

    ag.complete();
    return ag;
  }

  ArchGraph construct_vcol() {
    /*
     * 1 -- 1 -- 2  P1 -- C -- P2
     * |         |  |          |
     * 4         2  C          C
     * |         |  |          |
     * 4 -- 3 -- 3  P2 -- C -- P1
     */
    ArchGraph ag;

    auto p1 = ag.new_processor_type("P1");
    auto p2 = ag.new_processor_type("P2");
    auto c = ag.new_channel_type("C");

    auto pe1 = ag.add_processor(p1);
    auto pe2 = ag.add_processor(p2);
    auto pe3 = ag.add_processor(p1);
    auto pe4 = ag.add_processor(p2);

    ag.add_channel(pe1, pe2, c);
    ag.add_channel(pe2, pe3, c);
    ag.add_channel(pe3, pe4, c);
    ag.add_channel(pe4, pe1, c);

    ag.complete();
    return ag;
  }

  ArchGraph construct_ecol() {
    /*
     * 1 -- 1 -- 2  P -- C1 -- P
     * |         |  |          |
     * 4         2  C2         C2
     * |         |  |          |
     * 4 -- 3 -- 3  P -- C1 -- P
     */
    ArchGraph ag;

    auto p = ag.new_processor_type("P");
    auto c1 = ag.new_channel_type("C1");
    auto c2 = ag.new_channel_type("C2");

    auto pe1 = ag.add_processor(p);
    auto pe2 = ag.add_processor(p);
    auto pe3 = ag.add_processor(p);
    auto pe4 = ag.add_processor(p);

    ag.add_channel(pe1, pe2, c1);
    ag.add_channel(pe2, pe3, c2);
    ag.add_channel(pe3, pe4, c1);
    ag.add_channel(pe4, pe1, c2);

    ag.complete();
    return ag;
  }

  ArchGraph construct_tcol() {
    /*
     * 1 -- 1 -- 2  P1 -- C1 -- P2
     * |         |  |           |
     * 4         2  C2          C2
     * |         |  |           |
     * 4 -- 3 -- 3  P2 -- C1 -- P1
     */
    ArchGraph ag;

    auto p1 = ag.new_processor_type("P1");
    auto p2 = ag.new_processor_type("P2");
    auto c1 = ag.new_channel_type("C1");
    auto c2 = ag.new_channel_type("C2");

    auto pe1 = ag.add_processor(p1);
    auto pe2 = ag.add_processor(p2);
    auto pe3 = ag.add_processor(p1);
    auto pe4 = ag.add_processor(p2);

    ag.add_channel(pe1, pe2, c1);
    ag.add_channel(pe2, pe3, c2);
    ag.add_channel(pe3, pe4, c1);
    ag.add_channel(pe4, pe1, c2);

    ag.complete();
    return ag;
  }
};

class ArchGraphTest : public ArchGraphTestBase<testing::Test>{};

TEST_F(ArchGraphTest, CanObtainAutomorphisms)
{
  EXPECT_TRUE(perm_group_equal({
      {{1, 2, 3, 4}}, {{1, 3}, {2, 4}}, {{1, 4, 3, 2}}, {{1, 4}, {2, 3}},
      {{1, 2}, {3, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_nocol.automorphisms()))
    << "Automorphisms of uncolored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_vcol.automorphisms()))
    << "Automorphisms of processor colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}, {{1, 4}, {2, 3}}, {{1, 2}, {3, 4}}
    },
    ag_ecol.automorphisms()))
    << "Automorphisms of channel colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}
    }, ag_tcol.automorphisms()))
    << "Automorphisms of totally colored architecture graph correct.";
}

TEST_F(ArchGraphTest, CanLoadFromLua)
{
  ArchGraph ag = ArchGraph::fromlua(resource_path("mcsoc.lua"));

  EXPECT_EQ(8u, ag.num_processors())
    << "Loaded architecture graph has correct number of processors.";

  EXPECT_EQ(64u, ag.num_channels())
    << "Loaded architecture graph has correct number of channels.";
}

TEST_F(ArchGraphTest, CanProduceDotFile)
{
  ArchGraph ag = ArchGraph::fromlua(resource_path("mcsoc.lua"));

  ag.todot(resource_path("mcsoc.dot"));
}

class ArchGraphMappingVariantTest :
  public ArchGraphTestBase<
    testing::TestWithParam<ArchGraph::MappingVariant>> {};

TEST_P(ArchGraphMappingVariantTest, CanTestMappingEquivalence)
{
  std::vector<ArchGraphSystem const *> const arch_graphs {
    &ag_nocol, &ag_vcol, &ag_ecol, &ag_tcol
  };

  std::vector<std::vector<orbit>> const expected_orbits {
    {
      {{0, 0}, {1, 1}, {2, 2}, {3, 3}},
      {{0, 1}, {0, 3}, {1, 0}, {1, 2}, {2, 1}, {2, 3}, {3, 0}, {3, 2}},
      {{0, 2}, {1, 3}, {2, 0}, {3, 1}}
    },
    {
      {{0, 0}, {2, 2}},
      {{0, 1}, {0, 3}, {2, 1}, {2, 3}},
      {{0, 2}, {2, 0}},
      {{1, 0}, {1, 2}, {3, 0}, {3, 2}},
      {{1, 1}, {3, 3}},
      {{1, 3}, {3, 1}}
    },
    {
      {{0, 0}, {1, 1}, {2, 2}, {3, 3}},
      {{0, 1}, {1, 0}, {2, 3}, {3, 2}},
      {{0, 2}, {1, 3}, {2, 0}, {3, 1}},
      {{0, 3}, {1, 2}, {2, 1}, {3, 0}}
    },
    {
      {{0, 0}, {2, 2}},
      {{0, 1}, {2, 3}},
      {{0, 2}, {2, 0}},
      {{0, 3}, {2, 1}},
      {{1, 0}, {3, 2}},
      {{1, 1}, {3, 3}},
      {{1, 2}, {3, 0}},
      {{1, 3}, {3, 1}}
    },
  };

  for (auto i = 0u; i < arch_graphs.size(); ++i) {
    expect_mapping_generates_orbits(
      arch_graphs[i], expected_orbits[i], GetParam());
  }
}

INSTANTIATE_TEST_CASE_P(
  ArchGraphMappingVariants, ArchGraphMappingVariantTest,
  testing::Values(ArchGraphSystem::MAP_BRUTEFORCE,
                  ArchGraphSystem::MAP_APPROX));

template<typename T>
class ArchGraphClusterTestBase : public T
{
protected:
  void SetUp() {
    cluster_minimal = construct_minimal();
  }

  ArchGraphCluster cluster_minimal;

private:
  ArchGraphCluster construct_minimal() {
    /*
     * 1 -- 1 -- 2
     *
     * P -- C -- P
     */
    ArchGraph ag;

    auto p = ag.new_processor_type("P");
    auto c = ag.new_channel_type("C");

    auto pe1 = ag.add_processor(p);
    auto pe2 = ag.add_processor(p);

    ag.add_channel(pe1, pe2, c);

    /*
     * 1 -- 1 -- 2     3 -- 2 -- 4
     * |    |    |     |    |    |
     * ===========================
     *
     * P -- C -- P     P -- C -- P
     * |    |    |     |    |    |
     * ===========================
     */
    ArchGraphCluster cluster;

    cluster.add_subsystem(std::shared_ptr<ArchGraph>(new ArchGraph(ag)));
    cluster.add_subsystem(std::shared_ptr<ArchGraph>(new ArchGraph(ag)));

    cluster.complete();
    return cluster;
  }
};

class ArchGraphClusterTest : public ArchGraphClusterTestBase<testing::Test>{};

TEST_F(ArchGraphClusterTest, CanDetermineNumberOfProcessors)
{
  EXPECT_EQ(4u, cluster_minimal.num_processors())
    << "Number of processors in architecture graph cluster determined correctly.";
}

TEST_F(ArchGraphClusterTest, CanDetermineNumberOfChannels)
{
  EXPECT_EQ(2u, cluster_minimal.num_channels())
    << "Number of channels in architecture graph cluster determined correctly.";
}

class ArchGraphClusterMappingVariantTest :
  public ArchGraphClusterTestBase<
    testing::TestWithParam<ArchGraphSystem::MappingVariant>> {};

TEST_P(ArchGraphClusterMappingVariantTest, CanTestMappingEquivalence)
{
  std::vector<ArchGraphSystem const *> const clusters {&cluster_minimal};

  std::vector<std::vector<orbit>> const expected_orbits {
    {
      {{0, 0}, {1, 1}},
      {{0, 1}, {1, 0}},
      {{0, 2}, {0, 3}, {1, 2}, {1, 3}},
      {{2, 0}, {2, 1}, {3, 0}, {3, 1}},
      {{2, 2}, {3, 3}},
      {{2, 3}, {3, 2}}
    }
  };

  for (auto i = 0u; i < clusters.size(); ++i) {
    expect_mapping_generates_orbits(
      clusters[i], expected_orbits[i], GetParam());
  }
}

INSTANTIATE_TEST_CASE_P(
  ArchGraphClusterMappingVariants, ArchGraphClusterMappingVariantTest,
  testing::Values(ArchGraphSystem::MAP_BRUTEFORCE,
                  ArchGraphSystem::MAP_APPROX));
