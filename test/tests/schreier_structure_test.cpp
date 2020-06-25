#include <algorithm>
#include <memory>
#include <vector>

#include "gmock/gmock.h"

#include "explicit_transversals.hpp"
#include "orbits.hpp"
#include "perm.hpp"
#include "perm_set.hpp"
#include "schreier_tree.hpp"

#include "test_main.cpp"

using namespace mpsym;
using namespace mpsym::internal;

using testing::UnorderedElementsAreArray;

template <typename T>
class SchreierStructureTest : public testing::Test {};

using SchreierStructureTypes = ::testing::Types<ExplicitTransversals,
                                                SchreierTree>;

TYPED_TEST_SUITE(SchreierStructureTest, SchreierStructureTypes,);

TYPED_TEST(SchreierStructureTest, CanConstructSchreierStructures)
{
  unsigned n = 8;

  PermSet generators {
    Perm(n, {{1, 2, 3}}),
    Perm(n, {{1, 3}}),
    Perm(n, {{4, 6, 5}}),
    Perm(n, {{5, 6}, {7, 8}})
  };

  generators.insert_inverses();

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

  for (unsigned i = 0u; i < n; ++i) {
    unsigned root = i + 1u;

    auto schreier_structure(std::make_shared<TypeParam>(n, root, generators));

    Orbit::generate(root, generators, schreier_structure);

    EXPECT_EQ(root, schreier_structure->root())
      << "Root correct";

    EXPECT_THAT(expected_orbits[i],
                UnorderedElementsAreArray(schreier_structure->nodes()))
      << "Node (orbit) correct "
      << "(root is " << root << ").";

    for (unsigned x = 1u; x < n; ++x) {
      auto it(std::find(expected_orbits[i].begin(), expected_orbits[i].end(), x));
      bool contained = it != expected_orbits[i].end();

      EXPECT_EQ(contained, schreier_structure->contains(x))
        << "Can identify contained elements "
        << "(root is " << root << ", element is " << x << ").";
    }

    auto labels(schreier_structure->labels());

    std::vector<Perm> labels_vect(labels.begin(), labels.end());
    std::vector<Perm> gen_vect(generators.begin(), generators.end());

    EXPECT_THAT(labels_vect, UnorderedElementsAreArray(gen_vect))
      << "Edge labels correct "
      << "(root is " << root << ").";

    for (unsigned j = 0u; j < expected_orbits[i].size(); ++j) {
      unsigned origin = expected_orbits[i][j];

      Perm transv(schreier_structure->transversal(origin));

      EXPECT_EQ(origin, transv[root])
        << "Transversal " << transv << " correct "
        << "(root is " << root << ", origin is " << origin << ").";
    }
  }
}
