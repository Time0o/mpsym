#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>

#include "arch_graph.h"
#include "gmock/gmock.h"
#include "perm.h"
#include "perm_group.h"
#include "test_utility.h"

#include "test_main.cc"

using cgtl::ArchGraph;
using cgtl::Perm;
using cgtl::PermGroup;
using cgtl::TaskMapping;

using testing::UnorderedElementsAreArray;

class ArchGraphTest : public testing::Test
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

TEST_F(ArchGraphTest, CanTestMappingEquivalence)
{
  std::vector<ArchGraph> arch_graphs {ag_nocol, ag_vcol, ag_ecol, ag_tcol};

  typedef std::vector<std::vector<unsigned>> orbit;
  std::vector<std::vector<orbit>> expected_orbits {
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
    ArchGraph ag = arch_graphs[i];

    std::vector<TaskMapping> task_mappings;
    for (auto i = 0u; i < ag.num_processors(); ++i) {
      for (auto j = 0u; j < ag.num_processors(); ++j)
        task_mappings.push_back(ag.mapping({i, j}));
    }

    for (auto const &tm1 : task_mappings) {
      std::vector<unsigned> mapping1(tm1.mapping());

      std::stringstream ss; ss << "{ ";
      for (auto i = 0u; i < mapping1.size(); ++i)
        ss << mapping1[i] << (i == mapping1.size() - 1u ? " }" : ", ");

      std::vector<std::vector<unsigned>> equivalent_assignments;
      for (auto const &tm2 : task_mappings) {
        if (tm1.equivalence_class() == tm2.equivalence_class())
          equivalent_assignments.push_back(tm2.mapping());
      }

      bool found_orbit = false;
      for (auto const &orbit : expected_orbits[i]) {
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

  ag.todot("resources/mcsoc.dot");
}
