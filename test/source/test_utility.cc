#include <sstream>
#include <vector>

#include "gmock/gmock.h"
#include "perm.h"

template <typename P>
static ::testing::AssertionResult _perm_equal(
  std::vector<unsigned> const &expected, P const &perm)
{
  bool success = true;

  if (perm.degree() != expected.size()) {
    return ::testing::AssertionFailure()
      << "Permutation has incorrect degree (expected " << expected.size()
      << " but got " << perm.degree();
  }

  std::stringstream err;
  err << "Permutation differs:\n";

  for (unsigned i = 0u; i < perm.degree(); ++i) {
    if (perm[i + 1u] != expected[i]) {
      success = false;
      err << "@ index " << i + 1u << ":"
          << " expected " << expected[i]
          << " but got " << perm[i + 1u] << '\n';
    }
  }

  if (!success)
    return ::testing::AssertionFailure() << err.str();
  else
    return ::testing::AssertionSuccess();
}

::testing::AssertionResult perm_equal(std::vector<unsigned> const &expected,
  cgtl::Perm const &perm)
{
  return _perm_equal<cgtl::Perm>(expected, perm);
}

::testing::AssertionResult perm_word_equal(std::vector<unsigned> const &expected,
  cgtl::PermWord const &pw)
{
  return _perm_equal<cgtl::PermWord>(expected, pw);
}
