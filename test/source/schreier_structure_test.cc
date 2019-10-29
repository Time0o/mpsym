#include <memory>
#include <string>
#include <vector>

#include "gmock/gmock.h"

#include "schreier_sims.h"
#include "schreier_structure.h"

#include "test_main.cc"

using cgtl::ExplicitTransversals;
using cgtl::Perm;
using cgtl::SchreierTree;

using testing::UnorderedElementsAreArray;

template <typename T>
class SchreierStructureTest : public testing::Test {};

using SchreierStructureTypes = ::testing::Types<ExplicitTransversals,
                                                SchreierTree>;

TYPED_TEST_CASE(SchreierStructureTest, SchreierStructureTypes);

TYPED_TEST(SchreierStructureTest, CanConstructSchreierStructures)
{
  std::string display_type(typeid(TypeParam).name());

  unsigned n = 8;

  std::vector<Perm> const generators {
    Perm(n, {{1, 2, 3}}),
    Perm(n, {{1, 3}}),
    Perm(n, {{4, 6, 5}}),
    Perm(n, {{5, 6}, {7, 8}})
  };

  std::vector<unsigned> expected_orbits[] = {
    {1, 2, 3},
    {1, 2, 3},
    {1, 2, 3},
    {4, 5, 6},
    {4, 5, 6},
    {4, 5, 6},
    {7, 8},
    {7, 8}
  };

  std::vector<Perm> expected_transversals[] = {
    {
      Perm(n),
      Perm(n, {{1, 2, 3}}),
      Perm(n, {{1, 3}})
    },
    {
      Perm(n, {{1, 3, 2}}),
      Perm(n),
      Perm(n, {{1, 2, 3}})
    },
    {
      Perm(n, {{1, 2, 3}}),
      Perm(n, {{1, 3, 2}}),
      Perm(n)
    },
    {
      Perm(n),
      Perm(n, {{4, 5, 6}}),
      Perm(n, {{4, 6, 5}})
    },
    {
      Perm(n, {{4, 6, 5}}),
      Perm(n),
      Perm(n, {{5, 6}, {7, 8}})
    },
    {
      Perm(n, {{4, 5, 6}}),
      Perm(n, {{4, 6, 5}}),
      Perm(n)
    },
    {
      Perm(n),
      Perm(n, {{5, 6}, {7, 8}})
    },
    {
      Perm(n, {{5, 6}, {7, 8}}),
      Perm(n)
    }
  };

  for (unsigned i = 0u; i < n; ++i) {
    auto schreier_structure(std::make_shared<TypeParam>(n));

// TODO
//    orbit(i + 1u, generators, schreier_structure);
//
//    EXPECT_EQ(i + 1u, schreier_structure->root())
//      << "Root correct "
//      << "(type is " << display_type << ").";
//
//    EXPECT_THAT(expected_orbits[i],
//                UnorderedElementsAreArray(schreier_structure->nodes()))
//      << "Node (orbit) correct "
//      << "(root is " << i + 1u
//      << ", type is " << display_type << ").";
//
//    for (unsigned x = 1u; x < n; ++x) {
//      auto it(std::find(expected_orbits[i].begin(),
//                        expected_orbits[i].end(), x));
//
//      bool contained = it != expected_orbits[i].end();
//
//      EXPECT_EQ(contained, schreier_structure->contains(x))
//        << "Can identify contained elements "
//        << "(root is " << i + 1u
//        << ", element is " << x
//        << ", type is " << display_type << ").";
//    }
//
//    EXPECT_THAT(schreier_structure->labels(),
//                UnorderedElementsAreArray(generators))
//      << "Edge labels correct "
//      << "(root is " << i + 1u
//      << ", type is " << display_type << ").";
//
//    for (unsigned j = 0u; j < expected_orbits[i].size(); ++j) {
//      unsigned origin = expected_orbits[i][j];
//
//      EXPECT_EQ(expected_transversals[i][j],
//                schreier_structure->transversal(origin))
//        << "Transversal correct "
//        << "(root is " << i + 1u
//        << ", origin is " << origin
//        << ", type is " << display_type << ").";
//    }
  }
}
