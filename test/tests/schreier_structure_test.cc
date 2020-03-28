#include <memory>
#include <vector>

#include "gmock/gmock.h"

#include "explicit_transversals.h"
#include "orbits.h"
#include "perm.h"
#include "perm_set.h"
#include "schreier_tree.h"

#include "test_main.cc"

using mpsym::ExplicitTransversals;
using mpsym::Orbit;
using mpsym::Perm;
using mpsym::PermSet;
using mpsym::SchreierTree;

using testing::UnorderedElementsAreArray;

template <typename T>
class SchreierStructureTest : public testing::Test {};

using SchreierStructureTypes = ::testing::Types<ExplicitTransversals,
                                                SchreierTree>;

TYPED_TEST_CASE(SchreierStructureTest, SchreierStructureTypes);

TYPED_TEST(SchreierStructureTest, CanConstructSchreierStructures)
{
  unsigned n = 8;

  PermSet const generators {
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

  PermSet expected_transversals[] = {
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
    auto schreier_structure(std::make_shared<TypeParam>(n, i + 1u, generators));

    Orbit::generate(i + 1u, generators, schreier_structure);

    EXPECT_EQ(i + 1u, schreier_structure->root())
      << "Root correct";

    EXPECT_THAT(expected_orbits[i],
                UnorderedElementsAreArray(schreier_structure->nodes()))
      << "Node (orbit) correct "
      << "(root is " << i + 1u << ").";

    for (unsigned x = 1u; x < n; ++x) {
      auto it(std::find(expected_orbits[i].begin(), expected_orbits[i].end(), x));
      bool contained = it != expected_orbits[i].end();

      EXPECT_EQ(contained, schreier_structure->contains(x))
        << "Can identify contained elements "
        << "(root is " << i + 1u
        << ", element is " << x << ").";
    }

    auto labels(schreier_structure->labels());

    std::vector<Perm> labels_vect(labels.begin(), labels.end());
    std::vector<Perm> gen_vect(generators.begin(), generators.end());

    EXPECT_THAT(labels_vect, UnorderedElementsAreArray(gen_vect))
      << "Edge labels correct "
      << "(root is " << i + 1u << ").";

    for (unsigned j = 0u; j < expected_orbits[i].size(); ++j) {
      unsigned origin = expected_orbits[i][j];

      EXPECT_EQ(expected_transversals[i][j],
                schreier_structure->transversal(origin))
        << "Transversal correct "
        << "(root is " << i + 1u
        << ", origin is " << origin << ").";
    }
  }
}
