#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>

#include "gmock/gmock.h"

#include "arch_graph.h"
#include "partial_perm.h"
#include "perm.h"
#include "perm_group.h"
#include "test_utility.h"

#include "test_main.cc"

using cgtl::ArchGraphSystem;
using cgtl::ArchGraph;
using cgtl::ArchGraphCluster;
using cgtl::ArchUniformSuperGraph;
using cgtl::PartialPerm;
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

    ag.complete();
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

    ag.complete();
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

    ag.complete();
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

    ag.complete();
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

    ag.complete();
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

    ag.complete();
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

    ag.complete();
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

/*
TEST_F(ArchGraphTest, CanLoadFromLua)
{
  ArchGraph ag(ArchGraph::fromlua(resource_path("mcsoc.lua")));

  EXPECT_EQ(8u, ag.num_processors())
    << "Loaded architecture graph has correct number of processors.";

  EXPECT_EQ(64u, ag.num_channels())
    << "Loaded architecture graph has correct number of channels.";
}
*/

/*
TEST_F(ArchGraphTest, CanProduceDotFile)
{
  ArchGraph ag(ArchGraph::fromlua(resource_path("mcsoc.lua")));

  ag.todot(resource_path("mcsoc.dot"));
}
*/

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
  public ArchGraphTestBase<
    testing::TestWithParam<ArchGraph::MappingVariant>> {};

TEST_P(ArchGraphMappingVariantTest, CanTestMappingEquivalence)
{
  ArchGraph ag1(ag_nocol());
  ArchGraph ag2(ag_vcol());
  ArchGraph ag3(ag_ecol());
  ArchGraph ag4(ag_tcol());

  std::vector<ArchGraphSystem const *> const arch_graphs {
    &ag1, &ag2, &ag3, &ag4
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

TEST_F(ArchGraphClusterTest, CanObtainAutormorphisms)
{
  EXPECT_TRUE(perm_group_equal({
      Perm(4, {{1, 2}}),
      Perm(4, {{3, 4}}),
      Perm(4, {{1, 2}, {3, 4}})
    }, cluster_minimal.automorphisms()))
    << "Automorphisms of minimal architecture graph cluster correct.";
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

template<typename T>
class ArchUniformSuperGraphTestBase : public T
{
protected:
  void SetUp() {
    supergraph_minimal = construct_minimal();
  }

  ArchUniformSuperGraph supergraph_minimal;

private:
  ArchUniformSuperGraph construct_minimal() {
    // construct subsystem prototype
    ArchGraph proto;

    auto p = proto.new_processor_type();
    auto c = proto.new_channel_type();

    auto pe1 = proto.add_processor(p);
    auto pe2 = proto.add_processor(p);
    auto pe3 = proto.add_processor(p);

    proto.add_channel(pe1, pe2, c);
    proto.add_channel(pe2, pe3, c);
    proto.add_channel(pe3, pe1, c);

    // construct uniform supergraph
    ArchUniformSuperGraph supergraph(std::make_shared<ArchGraph>(proto));

    auto ssc = supergraph.new_subsystem_channel_type();

    auto ss1 = supergraph.add_subsystem();
    auto ss2 = supergraph.add_subsystem();
    auto ss3 = supergraph.add_subsystem();
    auto ss4 = supergraph.add_subsystem();

    supergraph.add_subsystem_channel(ss1, ss2, ssc);
    supergraph.add_subsystem_channel(ss2, ss3, ssc);
    supergraph.add_subsystem_channel(ss3, ss4, ssc);
    supergraph.add_subsystem_channel(ss4, ss1, ssc);

    supergraph.complete();

    return supergraph;
  }
};

class ArchUniformSuperGraphTest
  : public ArchUniformSuperGraphTestBase<testing::Test>{};

TEST_F(ArchUniformSuperGraphTest, CanDetermineNumberOfProcessors)
{
  EXPECT_EQ(12u, supergraph_minimal.num_processors())
    << "Number of processors in uniform architecture supergraph determined correctly.";
}

TEST_F(ArchUniformSuperGraphTest, CanDetermineNumberOfChannels)
{
  EXPECT_EQ(16u, supergraph_minimal.num_channels())
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

  EXPECT_EQ(expected_automorphisms, supergraph_minimal.automorphisms())
    << "Automorphisms of uniform architecture supergraph correct.";
}
