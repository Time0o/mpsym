#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <vector>

#include "gmock/gmock.h"

#include "arch_graph.h"
#include "arch_graph_cluster.h"
#include "arch_graph_system.h"
#include "arch_uniform_super_graph.h"
#include "partial_perm.h"
#include "perm.h"
#include "perm_group.h"
#include "task_allocation.h"
#include "test_utility.h"

#include "test_main.cc"

using cgtl::ArchGraph;
using cgtl::ArchGraphCluster;
using cgtl::ArchGraphSubsystem;
using cgtl::ArchGraphSystem;
using cgtl::ArchUniformSuperGraph;
using cgtl::PartialPerm;
using cgtl::Perm;
using cgtl::PermGroup;
using cgtl::TaskAllocation;

using testing::UnorderedElementsAreArray;

typedef std::vector<std::vector<unsigned>> orbit;

static void expect_mapping_generates_orbits(
  std::shared_ptr<ArchGraphSystem> const &ag,
  std::vector<orbit> expected_orbits,
  ArchGraphSystem::MappingMethod method)
{
  std::unordered_map<TaskAllocation, std::vector<TaskAllocation>> orbits;

  for (auto i = 1u; i <= ag->num_processors(); ++i) {
    for (auto j = 1u; j <= ag->num_processors(); ++j) {
      TaskAllocation allocation({i, j});
      TaskAllocation representative(ag->mapping(allocation, method));

      auto it = orbits.find(representative);
      if (it == orbits.end())
        orbits[representative] = {representative};

      if (allocation != representative)
        orbits[representative].push_back(allocation);
    }
  }

  for (auto const &orbit : orbits) {
    TaskAllocation representative(orbit.first);
    std::vector<TaskAllocation> actual_orbit(orbit.second);

    std::stringstream ss;
    ss << "{ " << representative[0];
    for (auto i = 0u; i < representative.size(); ++i)
      ss << representative[i] << ", " << representative[i];

    bool found_orbit = false;

    for (auto const &expected_orbit : expected_orbits) {
      auto it = std::find(expected_orbit.begin(),
                          expected_orbit.end(),
                          representative);

      if (it != expected_orbit.end()) {
        EXPECT_THAT(actual_orbit, UnorderedElementsAreArray(expected_orbit))
          << "Orbit with representative " << ss.str() << " correct.";

        found_orbit = true;
        break;
      }
    }

    EXPECT_TRUE(found_orbit) << "Orbit representative " << ss.str()
                             << " present in some orbit.";
  }
}

template<typename T>
class ArchGraphTestBase : public T
{
protected:
  ArchGraph ag_nocol() {
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

    return ag;
  }

  ArchGraph ag_vcol() {
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

    return ag;
  }

  ArchGraph ag_ecol() {
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

    return ag;
  }

  ArchGraph ag_tcol() {
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

    return ag;
  }

  ArchGraph ag_tri() {
    ArchGraph ag;
    auto p = ag.new_processor_type("P");
    auto c = ag.new_channel_type("C");

    auto pe1 = ag.add_processor(p);
    auto pe2 = ag.add_processor(p);
    auto pe3 = ag.add_processor(p);

    ag.add_channel(pe1, pe2, c);
    ag.add_channel(pe2, pe3, c);
    ag.add_channel(pe3, pe1, c);

    return ag;
  }

  ArchGraph ag_grid22() {
    /*
     * P1--P2
     * |   |
     * P3--P4
     */
    ArchGraph ag;

    auto p = ag.new_processor_type("P");
    auto c = ag.new_channel_type("C");

    auto pe1 = ag.add_processor(p);
    auto pe2 = ag.add_processor(p);
    auto pe3 = ag.add_processor(p);
    auto pe4 = ag.add_processor(p);

    ag.add_channel(pe1, pe2, c);
    ag.add_channel(pe1, pe3, c);
    ag.add_channel(pe2, pe4, c);
    ag.add_channel(pe3, pe4, c);

    return ag;
  }

  ArchGraph ag_grid33() {
    /*
     * P1--P2--P3
     * |   |   |
     * P4--P5--P6
     * |   |   |
     * P7--P8--P9
     */
    ArchGraph ag;

    auto p = ag.new_processor_type("P");
    auto c = ag.new_channel_type("C");

    auto pe1 = ag.add_processor(p);
    auto pe2 = ag.add_processor(p);
    auto pe3 = ag.add_processor(p);
    auto pe4 = ag.add_processor(p);
    auto pe5 = ag.add_processor(p);
    auto pe6 = ag.add_processor(p);
    auto pe7 = ag.add_processor(p);
    auto pe8 = ag.add_processor(p);
    auto pe9 = ag.add_processor(p);

    ag.add_channel(pe1, pe2, c);
    ag.add_channel(pe1, pe4, c);
    ag.add_channel(pe2, pe3, c);
    ag.add_channel(pe2, pe5, c);
    ag.add_channel(pe3, pe6, c);
    ag.add_channel(pe4, pe5, c);
    ag.add_channel(pe4, pe7, c);
    ag.add_channel(pe5, pe6, c);
    ag.add_channel(pe5, pe8, c);
    ag.add_channel(pe6, pe9, c);
    ag.add_channel(pe7, pe8, c);
    ag.add_channel(pe8, pe9, c);

    return ag;
  }
};

class ArchGraphTest : public ArchGraphTestBase<testing::Test>{};

TEST_F(ArchGraphTest, CanObtainAutomorphisms)
{
  EXPECT_TRUE(perm_group_equal({
      Perm(4, {{1, 2, 3, 4}}),
      Perm(4, {{1, 3}, {2, 4}}),
      Perm(4, {{1, 4, 3, 2}}),
      Perm(4, {{1, 4}, {2, 3}}),
      Perm(4, {{1, 2}, {3, 4}}),
      Perm(4, {{1, 3}}),
      Perm(4, {{2, 4}})
    }, ag_nocol().automorphisms()))
    << "Automorphisms of uncolored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      Perm(4, {{1, 3}, {2, 4}}),
      Perm(4, {{1, 3}}),
      Perm(4, {{2, 4}})
    }, ag_vcol().automorphisms()))
    << "Automorphisms of processor colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      Perm(4, {{1, 3}, {2, 4}}),
      Perm(4, {{1, 4}, {2, 3}}),
      Perm(4, {{1, 2}, {3, 4}})
    },
    ag_ecol().automorphisms()))
    << "Automorphisms of channel colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      Perm(4, {{1, 3}, {2, 4}})
    }, ag_tcol().automorphisms()))
    << "Automorphisms of totally colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      Perm(3, {{1, 2, 3}}),
      Perm(3, {{1, 2}}),
      Perm(3, {{1, 3, 2}}),
      Perm(3, {{1, 3}}),
      Perm(3, {{2, 3}})
    }, ag_tri().automorphisms()))
    << "Automorphisms of minimal triangular architecture graph correct.";
}

TEST_F(ArchGraphTest, DISABLED_CanObtainPartialAutomorphisms)
{
  auto ag(ag_grid22());
  auto partial_perm_inverse_semigroup(ag.partial_automorphisms());

  FAIL() << "TODO";
}

TEST(SpecialArchGraphTest, CanConstructFullyConnected)
{
  for (unsigned i = 1u; i < 5u; ++i) {
    EXPECT_EQ(PermGroup::symmetric(i),
              ArchGraph::fully_connected(i).automorphisms())
      << "Fully connected architecture graph with " << i
      << " processing elements has correct automorphism group.";
  }
}

TEST(SpecialArchGraphTest, CanConstructRegularMesh)
{
  // TODO: test non-quadratic meshes
  for (unsigned i = 1u; i < 5u; ++i) {
    EXPECT_EQ(PermGroup::dihedral(8),
              ArchGraph::regular_mesh(i, i).automorphisms())
      << "Regular mesh architecture graph with " << i * i
      << " processing elements has correct automorphism group.";
  }
}

class ArchGraphMappingVariantTest :
  public ArchGraphTestBase<testing::TestWithParam<ArchGraphSystem::MappingMethod>>
{};

TEST_P(ArchGraphMappingVariantTest, CanTestMappingEquivalence)
{
  std::vector<std::shared_ptr<ArchGraphSystem>> const arch_graphs {
    std::make_shared<ArchGraph>(ag_nocol()),
    std::make_shared<ArchGraph>(ag_vcol()),
    std::make_shared<ArchGraph>(ag_ecol()),
    std::make_shared<ArchGraph>(ag_tcol())
  };

  std::vector<std::vector<orbit>> const expected_orbits {
    {
      {{1, 1}, {2, 2}, {3, 3}, {4, 4}},
      {{1, 2}, {1, 4}, {2, 1}, {2, 3}, {3, 2}, {3, 4}, {4, 1}, {4, 3}},
      {{1, 3}, {2, 4}, {3, 1}, {4, 2}}
    },
    {
      {{1, 1}, {3, 3}},
      {{1, 2}, {1, 4}, {3, 2}, {3, 4}},
      {{1, 3}, {3, 1}},
      {{2, 1}, {2, 3}, {4, 1}, {4, 3}},
      {{2, 2}, {4, 4}},
      {{2, 4}, {4, 2}}
    },
    {
      {{1, 1}, {2, 2}, {3, 3}, {4, 4}},
      {{1, 2}, {2, 1}, {3, 4}, {4, 3}},
      {{1, 3}, {2, 4}, {3, 1}, {4, 2}},
      {{1, 4}, {2, 3}, {3, 2}, {4, 1}}
    },
    {
      {{1, 1}, {3, 3}},
      {{1, 2}, {3, 4}},
      {{1, 3}, {3, 1}},
      {{1, 4}, {3, 2}},
      {{2, 1}, {4, 3}},
      {{2, 2}, {4, 4}},
      {{2, 3}, {4, 1}},
      {{2, 4}, {4, 2}}
    },
  };

  for (auto i = 0u; i < arch_graphs.size(); ++i) {
    expect_mapping_generates_orbits(
      arch_graphs[i], expected_orbits[i], GetParam());
  }
}

INSTANTIATE_TEST_CASE_P(
  ArchGraphMappingVariants,
  ArchGraphMappingVariantTest,
  testing::Values(ArchGraphSystem::MappingMethod::ITERATE,
                  ArchGraphSystem::MappingMethod::LOCAL_SEARCH,
                  ArchGraphSystem::MappingMethod::ORBITS));

template<typename T>
class ArchGraphClusterTestBase : public T
{
protected:
  void SetUp() {
    construct_minimal();
  }

  std::shared_ptr<ArchGraphCluster> cluster_minimal;

private:
  void construct_minimal() {
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
    cluster_minimal = std::make_shared<ArchGraphCluster>();
    cluster_minimal->add_subsystem(ag);
    cluster_minimal->add_subsystem(ag);
  }
};

class ArchGraphClusterTest : public ArchGraphClusterTestBase<testing::Test>{};

TEST_F(ArchGraphClusterTest, CanDetermineNumberOfProcessors)
{
  EXPECT_EQ(4u, cluster_minimal->num_processors())
    << "Number of processors in architecture graph cluster determined correctly.";
}

TEST_F(ArchGraphClusterTest, CanDetermineNumberOfChannels)
{
  EXPECT_EQ(2u, cluster_minimal->num_channels())
    << "Number of channels in architecture graph cluster determined correctly.";
}

TEST_F(ArchGraphClusterTest, CanObtainAutormorphisms)
{
  EXPECT_TRUE(perm_group_equal({
      Perm(4, {{1, 2}}),
      Perm(4, {{3, 4}}),
      Perm(4, {{1, 2}, {3, 4}})
    }, cluster_minimal->automorphisms()))
    << "Automorphisms of minimal architecture graph cluster correct.";
}

class ArchGraphClusterMappingVariantTest :
  public ArchGraphClusterTestBase<testing::TestWithParam<ArchGraphSystem::MappingMethod>>
{};

TEST_P(ArchGraphClusterMappingVariantTest, CanTestMappingEquivalence)
{
  std::vector<std::shared_ptr<ArchGraphSystem>> const clusters {cluster_minimal};

  std::vector<std::vector<orbit>> const expected_orbits {
    {
      {{1, 1}, {2, 2}},
      {{1, 2}, {2, 1}},
      {{1, 3}, {1, 4}, {2, 3}, {2, 4}},
      {{3, 1}, {3, 2}, {4, 1}, {4, 2}},
      {{3, 3}, {4, 4}},
      {{3, 4}, {4, 3}}
    }
  };

  for (auto i = 0u; i < clusters.size(); ++i) {
    expect_mapping_generates_orbits(
      clusters[i], expected_orbits[i], GetParam());
  }
}

INSTANTIATE_TEST_CASE_P(
  ArchGraphClusterMappingVariants,
  ArchGraphClusterMappingVariantTest,
  testing::Values(ArchGraphSystem::MappingMethod::ITERATE,
                  ArchGraphSystem::MappingMethod::LOCAL_SEARCH,
                  ArchGraphSystem::MappingMethod::ORBITS));

template<typename T>
class ArchUniformSuperGraphTestBase : public T
{
protected:
  void SetUp() {
    construct_minimal();
  }

  std::shared_ptr<ArchUniformSuperGraph> supergraph_minimal;

private:
  void construct_minimal() {
    // construct subsystem prototype
    ArchGraph ag;

    auto p = ag.new_processor_type();
    auto c = ag.new_channel_type();

    auto pe1 = ag.add_processor(p);
    auto pe2 = ag.add_processor(p);
    auto pe3 = ag.add_processor(p);

    ag.add_channel(pe1, pe2, c);
    ag.add_channel(pe2, pe3, c);
    ag.add_channel(pe3, pe1, c);

    // construct uniform supergraph
    supergraph_minimal = std::make_shared<ArchUniformSuperGraph>(ag);

    auto ssc = supergraph_minimal->new_subsystem_channel_type();

    auto ss1 = supergraph_minimal->add_subsystem();
    auto ss2 = supergraph_minimal->add_subsystem();
    auto ss3 = supergraph_minimal->add_subsystem();
    auto ss4 = supergraph_minimal->add_subsystem();

    supergraph_minimal->add_subsystem_channel(ss1, ss2, ssc);
    supergraph_minimal->add_subsystem_channel(ss2, ss3, ssc);
    supergraph_minimal->add_subsystem_channel(ss3, ss4, ssc);
    supergraph_minimal->add_subsystem_channel(ss4, ss1, ssc);
  }
};

class ArchUniformSuperGraphTest
  : public ArchUniformSuperGraphTestBase<testing::Test>{};

TEST_F(ArchUniformSuperGraphTest, CanDetermineNumberOfProcessors)
{
  EXPECT_EQ(12u, supergraph_minimal->num_processors())
    << "Number of processors in uniform architecture supergraph determined correctly.";
}

TEST_F(ArchUniformSuperGraphTest, CanDetermineNumberOfChannels)
{
  EXPECT_EQ(16u, supergraph_minimal->num_channels())
    << "Number of channels in uniform architecture supergraph determined correctly.";
}

TEST_F(ArchUniformSuperGraphTest, CanObtainAutormorphisms)
{
  PermGroup expected_automorphisms(12,
    {
      Perm(12, {{1, 2}}),
      Perm(12, {{1, 4, 7, 10}, {2, 5, 8, 11}, {3, 6, 9, 12}}),
      Perm(12, {{10, 11}}),
      Perm(12, {{11, 12}}),
      Perm(12, {{2, 3}}),
      Perm(12, {{4, 10}, {5, 11}, {6, 12}}),
      Perm(12, {{4, 5}}),
      Perm(12, {{5, 6}}),
      Perm(12, {{7, 8}}),
      Perm(12, {{8, 9}})
    }
  );

  EXPECT_EQ(expected_automorphisms, supergraph_minimal->automorphisms())
    << "Automorphisms of uniform architecture supergraph correct.";
}
