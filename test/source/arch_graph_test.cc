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

class ArchGraphTest : public ::testing::Test
{
protected:
  void SetUp() {
    ag_nocol = construct_nocol();
    ag_vcol = construct_vcol();
    ag_ecol = construct_ecol();
    ag_tcol = construct_tcol();

    ag_nocol.todot("ag_nocol.dot");
    ag_vcol.todot("ag_vcol.dot");
    ag_ecol.todot("ag_ecol.dot");
    ag_tcol.todot("ag_tcol.dot");
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

    return ag;
  }
};

TEST_F(ArchGraphTest, CanObtainVertexAutomorphisms)
{
  EXPECT_TRUE(perm_group_equal({
      {{1, 2, 3, 4}}, {{1, 3}, {2, 4}}, {{1, 4, 3, 2}}, {{1, 4}, {2, 3}},
      {{1, 2}, {3, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_nocol.automorphisms(ArchGraph::AUTOM_PROCESSORS)))
    << "Processor automorphisms of uncolored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_vcol.automorphisms(ArchGraph::AUTOM_PROCESSORS)))
    << "Processor automorphisms of processor colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 2, 3, 4}}, {{1, 3}, {2, 4}}, {{1, 4, 3, 2}}, {{1, 4}, {2, 3}},
      {{1, 2}, {3, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_ecol.automorphisms(ArchGraph::AUTOM_PROCESSORS)))
    << "Processor automorphisms of channel colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_vcol.automorphisms(ArchGraph::AUTOM_PROCESSORS)))
    << "Processor automorphisms of totally colored architecture graph correct.";
}

TEST_F(ArchGraphTest, CanObtainEdgeAutomorphisms)
{
  EXPECT_TRUE(perm_group_equal({
      {{1, 2, 3, 4}}, {{1, 3}, {2, 4}}, {{1, 4, 3, 2}}, {{1, 4}, {2, 3}},
      {{1, 2}, {3, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_nocol.automorphisms(ArchGraph::AUTOM_CHANNELS)))
    << "Channel automorphisms of uncolored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 2, 3, 4}}, {{1, 3}, {2, 4}}, {{1, 4, 3, 2}}, {{1, 4}, {2, 3}},
      {{1, 2}, {3, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_vcol.automorphisms(ArchGraph::AUTOM_CHANNELS)))
    << "Channel automorphisms of processor colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}, {{1, 4}, {2, 3}}, {{1, 2}, {3, 4}}
    },
    ag_ecol.automorphisms(ArchGraph::AUTOM_CHANNELS)))
    << "Channel automorphisms of channel colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}, {{1, 4}, {2, 3}}, {{1, 2}, {3, 4}}
    }, ag_tcol.automorphisms(ArchGraph::AUTOM_CHANNELS)))
    << "Channel automorphisms of totally colored architecture graph correct.";
}

TEST_F(ArchGraphTest, CanObtainTotalAutomorphisms)
{
  EXPECT_TRUE(perm_group_equal({
      {{1, 2, 3, 4}}, {{1, 3}, {2, 4}}, {{1, 4, 3, 2}}, {{1, 4}, {2, 3}},
      {{1, 2}, {3, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_nocol.automorphisms(ArchGraph::AUTOM_TOTAL)))
    << "Total automorphisms of uncolored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}, {{1, 3}}, {{2, 4}}
    }, ag_vcol.automorphisms(ArchGraph::AUTOM_TOTAL)))
    << "Total automorphisms of processor colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}, {{1, 4}, {2, 3}}, {{1, 2}, {3, 4}}
    },
    ag_ecol.automorphisms(ArchGraph::AUTOM_TOTAL)))
    << "Total automorphisms of channel colored architecture graph correct.";

  EXPECT_TRUE(perm_group_equal({
      {{1, 3}, {2, 4}}
    }, ag_tcol.automorphisms(ArchGraph::AUTOM_TOTAL)))
    << "Channel automorphisms of totally colored architecture graph correct.";
}

TEST_F(ArchGraphTest, CanLoadFromLua)
{
  ArchGraph ag;

  ASSERT_TRUE(ag.fromlua("resources/mcsoc.lua"))
    << "Can load architecture graph from lua description";
}
