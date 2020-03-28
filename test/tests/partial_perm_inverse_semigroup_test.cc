#include "partial_perm.h"
#include "partial_perm_inverse_semigroup.h"

#include "test_main.cc"

using mpsym::PartialPerm;
using mpsym::PartialPermInverseSemigroup;

class PartialPermInverseSemigroupTest : public testing::Test
{
protected:
  void SetUp() {
    std::vector<PartialPerm> const generators {
      PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {4, 6, 8, 1, 5, 2, 7, 3, 9}),
      PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {5, 7, 9, 2, 4, 1, 6, 3, 8}),
      PartialPerm({2, 5, 6}, {5, 6, 2}),
      PartialPerm({1, 2, 3}, {3, 1, 2})
    };

    inverse_semigroup = PartialPermInverseSemigroup(generators);
  }

  PartialPermInverseSemigroup inverse_semigroup;

  std::vector<PartialPerm> const expected_elements {
    PartialPerm({}, {}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {1, 2, 3, 4, 5, 6, 7, 8, 9}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {1, 2, 3, 7, 6, 5, 4, 9, 8}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {2, 1, 3, 5, 4, 7, 6, 9, 8}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {2, 1, 3, 6, 7, 4, 5, 8, 9}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {4, 6, 8, 1, 5, 2, 7, 3, 9}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {4, 6, 8, 7, 2, 5, 1, 9, 3}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {5, 7, 9, 2, 4, 1, 6, 3, 8}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {5, 7, 9, 6, 1, 4, 2, 8, 3}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {6, 4, 8, 2, 7, 1, 5, 3, 9}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {6, 4, 8, 5, 1, 7, 2, 9, 3}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {7, 5, 9, 1, 6, 2, 4, 3, 8}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {7, 5, 9, 4, 2, 6, 1, 8, 3}),
    PartialPerm({1, 2, 3}, {1, 2, 3}),
    PartialPerm({1, 2, 3}, {1, 3, 2}),
    PartialPerm({1, 2, 3}, {2, 1, 3}),
    PartialPerm({1, 2, 3}, {2, 3, 1}),
    PartialPerm({1, 2, 3}, {3, 1, 2}),
    PartialPerm({1, 2, 3}, {3, 2, 1}),
    PartialPerm({1, 2, 3}, {4, 6, 8}),
    PartialPerm({1, 2, 3}, {4, 8, 6}),
    PartialPerm({1, 2, 3}, {5, 7, 9}),
    PartialPerm({1, 2, 3}, {5, 9, 7}),
    PartialPerm({1, 2, 3}, {6, 4, 8}),
    PartialPerm({1, 2, 3}, {6, 8, 4}),
    PartialPerm({1, 2, 3}, {7, 5, 9}),
    PartialPerm({1, 2, 3}, {7, 9, 5}),
    PartialPerm({1, 2, 3}, {8, 4, 6}),
    PartialPerm({1, 2, 3}, {8, 6, 4}),
    PartialPerm({1, 2, 3}, {9, 5, 7}),
    PartialPerm({1, 2, 3}, {9, 7, 5}),
    PartialPerm({1, 4, 7}, {1, 4, 7}),
    PartialPerm({1, 4, 7}, {1, 7, 4}),
    PartialPerm({1, 4, 7}, {2, 5, 6}),
    PartialPerm({1, 4, 7}, {2, 6, 5}),
    PartialPerm({1, 4, 7}, {4, 1, 7}),
    PartialPerm({1, 4, 7}, {4, 7, 1}),
    PartialPerm({1, 4, 7}, {5, 2, 6}),
    PartialPerm({1, 4, 7}, {5, 6, 2}),
    PartialPerm({1, 4, 7}, {6, 2, 5}),
    PartialPerm({1, 4, 7}, {6, 5, 2}),
    PartialPerm({1, 4, 7}, {7, 1, 4}),
    PartialPerm({1, 4, 7}, {7, 4, 1}),
    PartialPerm({1}, {1}),
    PartialPerm({1}, {2}),
    PartialPerm({1}, {3}),
    PartialPerm({1}, {4}),
    PartialPerm({1}, {5}),
    PartialPerm({1}, {6}),
    PartialPerm({1}, {7}),
    PartialPerm({1}, {8}),
    PartialPerm({1}, {9}),
    PartialPerm({2, 5, 6}, {1, 4, 7}),
    PartialPerm({2, 5, 6}, {1, 7, 4}),
    PartialPerm({2, 5, 6}, {2, 5, 6}),
    PartialPerm({2, 5, 6}, {2, 6, 5}),
    PartialPerm({2, 5, 6}, {4, 1, 7}),
    PartialPerm({2, 5, 6}, {4, 7, 1}),
    PartialPerm({2, 5, 6}, {5, 2, 6}),
    PartialPerm({2, 5, 6}, {5, 6, 2}),
    PartialPerm({2, 5, 6}, {6, 2, 5}),
    PartialPerm({2, 5, 6}, {6, 5, 2}),
    PartialPerm({2, 5, 6}, {7, 1, 4}),
    PartialPerm({2, 5, 6}, {7, 4, 1}),
    PartialPerm({2}, {1}),
    PartialPerm({2}, {2}),
    PartialPerm({2}, {3}),
    PartialPerm({2}, {4}),
    PartialPerm({2}, {5}),
    PartialPerm({2}, {6}),
    PartialPerm({2}, {7}),
    PartialPerm({2}, {8}),
    PartialPerm({2}, {9}),
    PartialPerm({3}, {1}),
    PartialPerm({3}, {2}),
    PartialPerm({3}, {3}),
    PartialPerm({3}, {4}),
    PartialPerm({3}, {5}),
    PartialPerm({3}, {6}),
    PartialPerm({3}, {7}),
    PartialPerm({3}, {8}),
    PartialPerm({3}, {9}),
    PartialPerm({4, 6, 8}, {1, 2, 3}),
    PartialPerm({4, 6, 8}, {1, 3, 2}),
    PartialPerm({4, 6, 8}, {2, 1, 3}),
    PartialPerm({4, 6, 8}, {2, 3, 1}),
    PartialPerm({4, 6, 8}, {3, 1, 2}),
    PartialPerm({4, 6, 8}, {3, 2, 1}),
    PartialPerm({4, 6, 8}, {4, 6, 8}),
    PartialPerm({4, 6, 8}, {4, 8, 6}),
    PartialPerm({4, 6, 8}, {5, 7, 9}),
    PartialPerm({4, 6, 8}, {5, 9, 7}),
    PartialPerm({4, 6, 8}, {6, 4, 8}),
    PartialPerm({4, 6, 8}, {6, 8, 4}),
    PartialPerm({4, 6, 8}, {7, 5, 9}),
    PartialPerm({4, 6, 8}, {7, 9, 5}),
    PartialPerm({4, 6, 8}, {8, 4, 6}),
    PartialPerm({4, 6, 8}, {8, 6, 4}),
    PartialPerm({4, 6, 8}, {9, 5, 7}),
    PartialPerm({4, 6, 8}, {9, 7, 5}),
    PartialPerm({4}, {1}),
    PartialPerm({4}, {2}),
    PartialPerm({4}, {3}),
    PartialPerm({4}, {4}),
    PartialPerm({4}, {5}),
    PartialPerm({4}, {6}),
    PartialPerm({4}, {7}),
    PartialPerm({4}, {8}),
    PartialPerm({4}, {9}),
    PartialPerm({5, 7, 9}, {1, 2, 3}),
    PartialPerm({5, 7, 9}, {1, 3, 2}),
    PartialPerm({5, 7, 9}, {2, 1, 3}),
    PartialPerm({5, 7, 9}, {2, 3, 1}),
    PartialPerm({5, 7, 9}, {3, 1, 2}),
    PartialPerm({5, 7, 9}, {3, 2, 1}),
    PartialPerm({5, 7, 9}, {4, 6, 8}),
    PartialPerm({5, 7, 9}, {4, 8, 6}),
    PartialPerm({5, 7, 9}, {5, 7, 9}),
    PartialPerm({5, 7, 9}, {5, 9, 7}),
    PartialPerm({5, 7, 9}, {6, 4, 8}),
    PartialPerm({5, 7, 9}, {6, 8, 4}),
    PartialPerm({5, 7, 9}, {7, 5, 9}),
    PartialPerm({5, 7, 9}, {7, 9, 5}),
    PartialPerm({5, 7, 9}, {8, 4, 6}),
    PartialPerm({5, 7, 9}, {8, 6, 4}),
    PartialPerm({5, 7, 9}, {9, 5, 7}),
    PartialPerm({5, 7, 9}, {9, 7, 5}),
    PartialPerm({5}, {1}),
    PartialPerm({5}, {2}),
    PartialPerm({5}, {3}),
    PartialPerm({5}, {4}),
    PartialPerm({5}, {5}),
    PartialPerm({5}, {6}),
    PartialPerm({5}, {7}),
    PartialPerm({5}, {8}),
    PartialPerm({5}, {9}),
    PartialPerm({6}, {1}),
    PartialPerm({6}, {2}),
    PartialPerm({6}, {3}),
    PartialPerm({6}, {4}),
    PartialPerm({6}, {5}),
    PartialPerm({6}, {6}),
    PartialPerm({6}, {7}),
    PartialPerm({6}, {8}),
    PartialPerm({6}, {9}),
    PartialPerm({7}, {1}),
    PartialPerm({7}, {2}),
    PartialPerm({7}, {3}),
    PartialPerm({7}, {4}),
    PartialPerm({7}, {5}),
    PartialPerm({7}, {6}),
    PartialPerm({7}, {7}),
    PartialPerm({7}, {8}),
    PartialPerm({7}, {9}),
    PartialPerm({8}, {1}),
    PartialPerm({8}, {2}),
    PartialPerm({8}, {3}),
    PartialPerm({8}, {4}),
    PartialPerm({8}, {5}),
    PartialPerm({8}, {6}),
    PartialPerm({8}, {7}),
    PartialPerm({8}, {8}),
    PartialPerm({8}, {9}),
    PartialPerm({9}, {1}),
    PartialPerm({9}, {2}),
    PartialPerm({9}, {3}),
    PartialPerm({9}, {4}),
    PartialPerm({9}, {5}),
    PartialPerm({9}, {6}),
    PartialPerm({9}, {7}),
    PartialPerm({9}, {8}),
    PartialPerm({9}, {9})
  };

  std::vector<PartialPerm> const expected_non_elements {
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {2, 1, 9, 3, 7, 6, 5, 4, 8}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {2, 5, 9, 4, 7, 1, 6, 3, 8}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {2, 9, 6, 1, 3, 8, 4, 5, 7}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {4, 2, 1, 5, 8, 3, 7, 6, 9}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {4, 3, 8, 1, 9, 2, 7, 6, 5}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {6, 4, 9, 5, 1, 7, 8, 2, 3}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {6, 7, 5, 2, 3, 8, 4, 9, 1}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {7, 3, 5, 6, 1, 9, 8, 2, 4}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {7, 8, 6, 5, 9, 2, 1, 4, 3}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {7, 9, 6, 3, 4, 5, 1, 8, 2}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {8, 4, 2, 5, 6, 3, 1, 9, 7}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {8, 9, 7, 1, 3, 5, 2, 4, 6}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 8}, {3, 4, 8, 9, 7, 5, 1, 2}),
    PartialPerm({1, 2, 3, 4, 5, 6, 7, 9}, {1, 5, 6, 9, 2, 4, 3, 8}),
    PartialPerm({1, 2, 3, 4, 5, 6, 8, 9}, {7, 3, 2, 9, 8, 1, 4, 5}),
    PartialPerm({1, 2, 3, 4, 5, 6, 8, 9}, {8, 2, 1, 4, 5, 6, 9, 3}),
    PartialPerm({1, 2, 3, 4, 5, 7, 8}, {2, 8, 6, 5, 9, 4, 7}),
    PartialPerm({1, 2, 3, 4, 6, 7, 8, 9}, {1, 3, 8, 7, 9, 6, 2, 5}),
    PartialPerm({1, 2, 3, 4, 6, 7, 8, 9}, {1, 5, 8, 6, 2, 4, 7, 9}),
    PartialPerm({1, 2, 3, 4, 6, 7, 8, 9}, {6, 2, 8, 4, 1, 5, 3, 9}),
    PartialPerm({1, 2, 3, 4, 6, 7, 8}, {8, 1, 5, 3, 7, 2, 9}),
    PartialPerm({1, 2, 3, 4, 6, 7}, {1, 6, 7, 2, 4, 9}),
    PartialPerm({1, 2, 3, 4, 7, 8, 9}, {6, 3, 5, 2, 8, 9, 4}),
    PartialPerm({1, 2, 3, 4, 7, 8, 9}, {8, 1, 5, 6, 4, 2, 3}),
    PartialPerm({1, 2, 3, 4, 7, 9}, {2, 6, 8, 4, 7, 3}),
    PartialPerm({1, 2, 3, 4, 9}, {3, 2, 7, 6, 8}),
    PartialPerm({1, 2, 3, 5, 6, 7, 8, 9}, {8, 6, 7, 9, 4, 2, 5, 3}),
    PartialPerm({1, 2, 3, 6, 7, 8, 9}, {2, 4, 5, 6, 1, 8, 9}),
    PartialPerm({1, 2, 4, 5, 6, 7, 8, 9}, {4, 5, 6, 9, 1, 3, 2, 8}),
    PartialPerm({1, 2, 4, 5, 6, 7, 8, 9}, {4, 6, 9, 7, 2, 8, 3, 1}),
    PartialPerm({1, 2, 4, 5, 6, 7, 8, 9}, {6, 4, 1, 5, 9, 3, 8, 2}),
    PartialPerm({1, 2, 4, 5, 6, 7}, {5, 3, 2, 4, 8, 6}),
    PartialPerm({1, 2, 4, 5, 8, 9}, {8, 9, 6, 1, 7, 3}),
    PartialPerm({1, 2, 4, 6, 7, 8, 9}, {5, 8, 2, 6, 4, 9, 1}),
    PartialPerm({1, 2, 4, 6, 7, 9}, {1, 6, 8, 9, 3, 4}),
    PartialPerm({1, 2, 4, 6}, {2, 9, 1, 6}),
    PartialPerm({1, 2, 4, 7, 9}, {2, 3, 6, 1, 4}),
    PartialPerm({1, 2, 5, 6, 9}, {2, 4, 8, 5, 7}),
    PartialPerm({1, 2, 5, 9}, {5, 3, 4, 9}),
    PartialPerm({1, 2, 6, 7, 8, 9}, {7, 6, 1, 4, 5, 3}),
    PartialPerm({1, 2, 7, 8, 9}, {4, 1, 2, 9, 8}),
    PartialPerm({1, 3, 4, 5, 6, 7, 8, 9}, {2, 7, 8, 5, 1, 3, 4, 6}),
    PartialPerm({1, 3, 4, 5, 6, 7, 8}, {2, 5, 7, 3, 9, 6, 1}),
    PartialPerm({1, 3, 4, 5, 7, 8, 9}, {1, 5, 6, 8, 7, 4, 9}),
    PartialPerm({1, 3, 4, 5, 7, 8, 9}, {1, 6, 9, 3, 2, 5, 8}),
    PartialPerm({1, 3, 4, 6, 9}, {6, 5, 3, 4, 7}),
    PartialPerm({1, 3, 5}, {2, 1, 4}),
    PartialPerm({1, 3, 6, 7, 9}, {8, 2, 9, 7, 4}),
    PartialPerm({1, 3, 8}, {9, 8, 6}),
    PartialPerm({1, 4, 5, 6, 7, 8, 9}, {8, 7, 9, 2, 1, 3, 4}),
    PartialPerm({1, 4, 5, 6, 7}, {6, 1, 3, 8, 9}),
    PartialPerm({1, 4, 5, 6, 8, 9}, {5, 9, 8, 1, 2, 6}),
    PartialPerm({1, 4, 5, 6, 8}, {1, 4, 5, 8, 7}),
    PartialPerm({1, 4, 6, 7, 9}, {3, 4, 6, 8, 7}),
    PartialPerm({1, 5, 7}, {9, 4, 3}),
    PartialPerm({1, 5, 8, 9}, {4, 3, 5, 6}),
    PartialPerm({1, 5, 9}, {5, 7, 9}),
    PartialPerm({1, 5}, {3, 7}),
    PartialPerm({1, 6, 7, 8, 9}, {5, 7, 9, 1, 8}),
    PartialPerm({1, 6}, {5, 1}),
    PartialPerm({1, 7, 8, 9}, {3, 6, 5, 7}),
    PartialPerm({2, 3, 4, 5, 6, 7, 8, 9}, {6, 5, 8, 1, 9, 2, 3, 4}),
    PartialPerm({2, 3, 4, 5, 6, 7, 8}, {6, 1, 2, 7, 8, 3, 9}),
    PartialPerm({2, 3, 4, 5, 6, 7, 9}, {9, 2, 8, 3, 5, 4, 7}),
    PartialPerm({2, 3, 4, 5, 6, 7}, {4, 3, 1, 2, 6, 7}),
    PartialPerm({2, 3, 4, 5, 6, 9}, {3, 2, 1, 5, 8, 4}),
    PartialPerm({2, 3, 4, 7, 8}, {2, 4, 1, 3, 8}),
    PartialPerm({2, 3, 5, 6, 9}, {5, 8, 3, 2, 7}),
    PartialPerm({2, 3, 6, 7}, {8, 4, 6, 2}),
    PartialPerm({2, 4, 5, 6, 7, 9}, {7, 9, 3, 8, 5, 4}),
    PartialPerm({2, 4, 5, 6, 8, 9}, {1, 8, 6, 4, 3, 5}),
    PartialPerm({2, 5, 6, 7}, {7, 4, 5, 9}),
    PartialPerm({2, 5, 6, 9}, {2, 6, 4, 8}),
    PartialPerm({2, 5, 6}, {5, 6, 7}),
    PartialPerm({2, 5, 7, 8, 9}, {3, 4, 6, 8, 7}),
    PartialPerm({2, 6, 7, 9}, {5, 7, 9, 3}),
    PartialPerm({2, 6, 7}, {8, 2, 1}),
    PartialPerm({2, 6}, {2, 4}),
    PartialPerm({2, 7, 8}, {3, 6, 5}),
    PartialPerm({2, 7, 9}, {5, 8, 7}),
    PartialPerm({2, 7}, {7, 4}),
    PartialPerm({2, 8}, {1, 6}),
    PartialPerm({3, 4, 5, 6, 7, 8, 9}, {3, 9, 5, 1, 2, 8, 7}),
    PartialPerm({3, 4, 5, 7, 8, 9}, {6, 9, 2, 7, 3, 5}),
    PartialPerm({3, 4, 7, 9}, {9, 5, 2, 7}),
    PartialPerm({3, 5, 6, 7, 8}, {3, 8, 6, 5, 9}),
    PartialPerm({3, 5, 6, 8, 9}, {3, 9, 8, 4, 2}),
    PartialPerm({3, 6, 8}, {3, 8, 9}),
    PartialPerm({3, 6}, {9, 2}),
    PartialPerm({3, 7, 8, 9}, {9, 8, 5, 7}),
    PartialPerm({3, 8}, {9, 5}),
    PartialPerm({3, 9}, {2, 6}),
    PartialPerm({4, 5, 7, 9}, {5, 9, 4, 6}),
    PartialPerm({4, 6}, {7, 1}),
    PartialPerm({4, 7, 9}, {8, 7, 4}),
    PartialPerm({5, 6, 7, 8, 9}, {8, 5, 6, 4, 2}),
    PartialPerm({5, 6, 8, 9}, {8, 1, 5, 9}),
    PartialPerm({5, 7}, {6, 1}),
    PartialPerm({7, 8, 9}, {3, 1, 9}),
    PartialPerm({8, 9}, {7, 1}),
    PartialPerm({10, 11, 12}, {1, 2, 3}),
    PartialPerm({1, 2, 3}, {10, 11, 12})
  };
};

TEST_F(PartialPermInverseSemigroupTest, CanTestMembership)
{
  for (PartialPerm const &pperm : expected_elements) {
    EXPECT_TRUE(inverse_semigroup.contains_element(pperm))
      << "Can recognize inverse semigroup element (" << pperm << ").";
  }

  for (PartialPerm const &pperm : expected_non_elements) {
    EXPECT_FALSE(inverse_semigroup.contains_element(pperm))
      << "Can recognize inverse semigroup non-element (" << pperm << ").";
  }
}

TEST_F(PartialPermInverseSemigroupTest, CanAdjoinGenerators)
{
  PartialPermInverseSemigroup inverse_semigroup;

  inverse_semigroup.adjoin_generators(
    {PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {4, 6, 8, 1, 5, 2, 7, 3, 9})});

  inverse_semigroup.adjoin_generators(
      {PartialPerm({1, 2, 3, 4, 5, 6, 7, 8, 9}, {5, 7, 9, 2, 4, 1, 6, 3, 8})});

  inverse_semigroup.adjoin_generators(
      {PartialPerm({2, 5, 6}, {5, 6, 2})});

  inverse_semigroup.adjoin_generators(
      {PartialPerm({1, 2, 3}, {3, 1, 2})});

  for (PartialPerm const &pperm : expected_elements) {
    EXPECT_TRUE(inverse_semigroup.contains_element(pperm))
      << "Stepwise constructed inverse semigroup contains correct elements ("
      << pperm << ").";
  }

  for (PartialPerm const &pperm : expected_non_elements) {
    EXPECT_FALSE(inverse_semigroup.contains_element(pperm))
      << "Stepwise constructed inverse semigroup does not contain additional members"
      << pperm << ").";
  }
}
