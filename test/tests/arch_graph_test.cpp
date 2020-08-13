#include <algorithm>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gmock/gmock.h"

#include "arch_graph.hpp"
#include "arch_graph_cluster.hpp"
#include "arch_graph_system.hpp"
#include "arch_uniform_super_graph.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "task_mapping.hpp"
#include "test_utility.hpp"

#include "test_main.cpp"

using namespace mpsym;
using namespace mpsym::internal;

using testing::UnorderedElementsAreArray;

typedef std::vector<std::vector<unsigned>> orbit;

static void expect_generates_orbits(
  std::shared_ptr<ArchGraphSystem> const &ag,
  std::vector<orbit> expected_orbits,
  ReprOptions::Method method)
{
  std::unordered_map<TaskMapping, std::vector<TaskMapping>> orbits;

  for (auto i = 1u; i <= ag->num_processors(); ++i) {
    for (auto j = 1u; j <= ag->num_processors(); ++j) {
      TaskMapping mapping({i, j});

      ReprOptions options;
      options.method = method;

      TaskMapping repr(ag->repr(mapping, nullptr, &options));

      if (orbits.find(repr) == orbits.end())
        orbits[repr] = {repr};

      if (mapping != repr)
        orbits[repr].push_back(mapping);
    }
  }

  for (auto const &orbit : orbits) {
    TaskMapping representative(orbit.first);
    std::vector<TaskMapping> actual_orbit(orbit.second);

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

class ArchGraphReprVariantTest :
  public ArchGraphTestBase<testing::TestWithParam<ReprOptions::Method>>
{};

TEST_P(ArchGraphReprVariantTest, CanTestReprEquivalence)
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

  for (auto i = 0u; i < arch_graphs.size(); ++i)
    expect_generates_orbits(arch_graphs[i], expected_orbits[i], GetParam());
}

INSTANTIATE_TEST_SUITE_P(
  ArchGraphReprVariants,
  ArchGraphReprVariantTest,
  testing::Values(ReprOptions::Method::ITERATE,
                  ReprOptions::Method::LOCAL_SEARCH,
                  ReprOptions::Method::ORBITS));

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
    auto ag(std::make_shared<ArchGraph>());

    auto p = ag->new_processor_type("P");
    auto c = ag->new_channel_type("C");

    auto pe1 = ag->add_processor(p);
    auto pe2 = ag->add_processor(p);

    ag->add_channel(pe1, pe2, c);

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

class ArchGraphClusterReprVariantTest :
  public ArchGraphClusterTestBase<testing::TestWithParam<ReprOptions::Method>>
{};

TEST_P(ArchGraphClusterReprVariantTest, CanTestReprEquivalence)
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

  for (auto i = 0u; i < clusters.size(); ++i)
    expect_generates_orbits(clusters[i], expected_orbits[i], GetParam());
}

INSTANTIATE_TEST_SUITE_P(
  ArchGraphClusterReprVariants,
  ArchGraphClusterReprVariantTest,
  testing::Values(ReprOptions::Method::ITERATE,
                  ReprOptions::Method::LOCAL_SEARCH,
                  ReprOptions::Method::ORBITS));

template<typename T>
class ArchUniformSuperGraphTestBase : public T
{
protected:
  void SetUp() {
    construct_minimal();
  }

  std::shared_ptr<ArchUniformSuperGraph> super_graph_minimal;

private:
  void construct_minimal() {
    // construct uniform super_graph
    auto super_graph(std::make_shared<ArchGraph>());

    auto sp = super_graph->new_processor_type();
    auto sc = super_graph->new_channel_type();

    auto spe1 = super_graph->add_processor(sp);
    auto spe2 = super_graph->add_processor(sp);
    auto spe3 = super_graph->add_processor(sp);
    auto spe4 = super_graph->add_processor(sp);

    super_graph->add_channel(spe1, spe2, sc);
    super_graph->add_channel(spe2, spe3, sc);
    super_graph->add_channel(spe3, spe4, sc);
    super_graph->add_channel(spe4, spe1, sc);

    // construct subsystem prototype
    auto proto(std::make_shared<ArchGraph>());

    auto p = proto->new_processor_type();
    auto c = proto->new_channel_type();

    auto pe1 = proto->add_processor(p);
    auto pe2 = proto->add_processor(p);
    auto pe3 = proto->add_processor(p);

    proto->add_channel(pe1, pe2, c);
    proto->add_channel(pe2, pe3, c);
    proto->add_channel(pe3, pe1, c);

    super_graph_minimal =
      std::make_shared<ArchUniformSuperGraph>(super_graph, proto);
  }
};

class ArchUniformSuperGraphTest
  : public ArchUniformSuperGraphTestBase<testing::Test>{};

TEST_F(ArchUniformSuperGraphTest, CanDetermineNumberOfProcessors)
{
  EXPECT_EQ(12u, super_graph_minimal->num_processors())
    << "Number of processors in uniform architecture super_graph determined correctly.";
}

TEST_F(ArchUniformSuperGraphTest, CanDetermineNumberOfChannels)
{
  EXPECT_EQ(16u, super_graph_minimal->num_channels())
    << "Number of channels in uniform architecture super_graph determined correctly.";
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

  EXPECT_EQ(expected_automorphisms, super_graph_minimal->automorphisms())
    << "Automorphisms of uniform architecture super_graph correct.";
}
