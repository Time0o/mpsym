#include <sstream>
#include <vector>

#include "gmock/gmock-matchers.h"
#include "gmock/gmock.h"
#include "perm.h"
#include "perm_group.h"

using cgtl::Perm;
using cgtl::PermWord;
using cgtl::PermGroup;

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
  Perm const &perm)
{
  return _perm_equal<Perm>(expected, perm);
}

::testing::AssertionResult perm_word_equal(std::vector<unsigned> const &expected,
  PermWord const &pw)
{
  return _perm_equal<PermWord>(expected, pw);
}

::testing::AssertionResult perm_group_equal(
  std::vector<std::vector<std::vector<unsigned>>> const &expected,
  PermGroup const &pg)
{
    unsigned degree = pg.degree();

    std::vector<Perm> expected_elements {Perm(degree)};
    for (auto const &p : expected) {
      expected_elements.push_back(Perm(degree, p));
    }

    std::vector<Perm> actual_elements;
    for (Perm const &p : pg)
      actual_elements.push_back(p);

    ::testing::Matcher<std::vector<Perm> const &> matcher(
      ::testing::UnorderedElementsAreArray(expected_elements));

    ::testing::StringMatchResultListener listener;
    if (matcher.MatchAndExplain(actual_elements, &listener))
      return ::testing::AssertionSuccess();

    std::stringstream ss;
    ss << "\nShould be: ";
    matcher.DescribeTo(&ss);

    ss << "\nBut is: { ";
    for (auto i = 0u; i < actual_elements.size(); ++i) {
      ss << actual_elements[i]
         << (i == actual_elements.size() - 1u ? "" : ", ");
    }
    ss << " },\n";
    ss << listener.str();

    std::string msg(ss.str());

    int const indent = 4;

    size_t i = 0u;
    for (;;) {
      if ((i = msg.find('\n', i)) == std::string::npos)
        break;

      msg.insert(i + 1, std::string(indent, ' '));
      i += indent + 1;
    }

    return ::testing::AssertionFailure() << msg;
}
