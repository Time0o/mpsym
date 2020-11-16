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
    Perm(n, {{0, 1, 2}}),
    Perm(n, {{0, 2}}),
    Perm(n, {{3, 5, 4}}),
    Perm(n, {{4, 5}, {6, 7}})
  };

  generators.insert_inverses();

  std::vector<unsigned> expected_orbits[] = {
    {0, 1, 2},
    {0, 1, 2},
    {0, 1, 2},
    {3, 4, 5},
    {3, 4, 5},
    {3, 4, 5},
    {6, 7},
    {6, 7}
  };

  for (unsigned root = 0u; root < n; ++root) {
    auto schreier_structure(std::make_shared<TypeParam>(n, root, generators));

    Orbit::generate(root, generators, schreier_structure);

    EXPECT_EQ(root, schreier_structure->root())
      << "Root correct";

    auto const &orbit(expected_orbits[root]);

    EXPECT_THAT(orbit, UnorderedElementsAreArray(schreier_structure->nodes()))
      << "Node (orbit) correct "
      << "(root is " << root << ").";

    for (unsigned x = 1u; x < n; ++x) {
      auto it(std::find(orbit.begin(), orbit.end(), x));
      bool contained = it != orbit.end();

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

    for (unsigned j = 0u; j < orbit.size(); ++j) {
      unsigned origin = orbit[j];

      Perm transv(schreier_structure->transversal(origin));

      EXPECT_EQ(origin, transv[root])
        << "Transversal " << transv << " correct "
        << "(root is " << root << ", origin is " << origin << ").";
    }
  }
}
