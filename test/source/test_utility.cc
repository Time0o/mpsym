#include <cassert>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

#include "gmock/gmock.h"
#include "gmock/gmock-matchers.h"

#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"
#include "test_utility.h"

using mpsym::Perm;
using mpsym::PermGroup;
using mpsym::PermSet;
using mpsym::PermWord;

template <typename P>
static testing::AssertionResult _perm_equal(
  std::vector<unsigned> const &expected, P const &perm)
{
  bool success = true;

  if (perm.degree() != expected.size()) {
    return testing::AssertionFailure()
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
    return testing::AssertionFailure() << err.str();
  else
    return testing::AssertionSuccess();
}

testing::AssertionResult perm_equal(std::vector<unsigned> const &expected,
  Perm const &perm)
{
  return _perm_equal<Perm>(expected, perm);
}

testing::AssertionResult perm_word_equal(std::vector<unsigned> const &expected,
  PermWord const &pw)
{
  return _perm_equal<PermWord>(expected, pw);
}

testing::AssertionResult perm_group_equal(PermGroup const &expected,
                                          PermGroup const &actual)
{
  PermSet expected_elements;
  for (Perm const &perm : expected)
    expected_elements.insert(perm);

  return perm_group_equal(expected_elements, actual);
}

testing::AssertionResult perm_group_equal(PermSet expected_elements,
                                          PermGroup const &actual)
{
  bool expected_has_id = false;
  for (Perm const &perm : expected_elements) {
    if (perm.id()) {
      expected_has_id = true;
      break;
    }
  }

  if (!expected_has_id)
    expected_elements.emplace(actual.degree());

  std::vector<Perm> actual_elements;
  for (Perm const &perm : actual)
    actual_elements.push_back(perm);

  testing::Matcher<std::vector<Perm> const &> matcher(
    testing::UnorderedElementsAreArray(expected_elements.begin(),
                                       expected_elements.end()));

  testing::StringMatchResultListener listener;
  if (matcher.MatchAndExplain(actual_elements, &listener))
    return testing::AssertionSuccess();

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

  return testing::AssertionFailure() << msg;
}

PermGroup verified_perm_group(VerifiedGroup group)
{
  struct PermutationGroupDescription {
    PermutationGroupDescription() {}

    PermutationGroupDescription(unsigned degree,
                                PermSet generators,
                                PermSet elements)
     : degree(degree),
       generators(generators),
       elements(elements) {}

    unsigned degree;
    PermSet generators;
    PermSet elements;
#ifndef NDEBUG
    bool verified = false;
#endif
  };

  static std::map<VerifiedGroup, PermutationGroupDescription> verified_groups {
    {S1,
      PermutationGroupDescription(1,
        {},
        {Perm(1)})
    },
    {S2,
      PermutationGroupDescription(2,
        {Perm(2, {{1, 2}})},
        {Perm(2),
         Perm(2, {{1, 2}})})
    },
    {S3,
      PermutationGroupDescription(3,
        {Perm(3, {{1, 2}}),
         Perm(3, {{1, 2, 3}})},
        {Perm(3),
         Perm(3, {{1, 2, 3}}),
         Perm(3, {{1, 2}}),
         Perm(3, {{1, 3, 2}}),
         Perm(3, {{1, 3}}),
         Perm(3, {{2, 3}})})
    },
    {S4,
      PermutationGroupDescription(4,
        {Perm(4, {{1, 2}}),
         Perm(4, {{1, 2, 3, 4}})},
        {Perm(4),
         Perm(4, {{1, 2, 3, 4}}),
         Perm(4, {{1, 2, 3}}),
         Perm(4, {{1, 2, 4, 3}}),
         Perm(4, {{1, 2, 4}}),
         Perm(4, {{1, 2}, {3,4}}),
         Perm(4, {{1, 2}}),
         Perm(4, {{1, 3, 2}}),
         Perm(4, {{1, 3, 2, 4}}),
         Perm(4, {{1, 3, 4, 2}}),
         Perm(4, {{1, 3, 4}}),
         Perm(4, {{1, 3}, {2, 4}}),
         Perm(4, {{1, 3}}),
         Perm(4, {{1, 4, 2, 3}}),
         Perm(4, {{1, 4, 2}}),
         Perm(4, {{1, 4, 3, 2}}),
         Perm(4, {{1, 4, 3}}),
         Perm(4, {{1, 4}, {2,3}}),
         Perm(4, {{1, 4}}),
         Perm(4, {{2, 3, 4}}),
         Perm(4, {{2, 3}}),
         Perm(4, {{2, 4, 3}}),
         Perm(4, {{2, 4}}),
         Perm(4, {{3, 4}})})
    },
    {S5,
      PermutationGroupDescription(5,
        {Perm(5, {{1, 2}}),
         Perm(5, {{1, 2, 3, 4, 5}})},
        {Perm(5),
         Perm(5, {{1, 2}, {3,4}}),
         Perm(5, {{1, 2}, {3,4,5}}),
         Perm(5, {{1, 2}, {3,5}}),
         Perm(5, {{1, 2}, {3,5,4}}),
         Perm(5, {{1, 2}, {4,5}}),
         Perm(5, {{1, 2}}),
         Perm(5, {{1, 2, 3}, {4, 5}}),
         Perm(5, {{1, 2, 3}}),
         Perm(5, {{1, 2, 3, 4}}),
         Perm(5, {{1, 2, 3, 4, 5}}),
         Perm(5, {{1, 2, 3, 5}}),
         Perm(5, {{1, 2, 3, 5, 4}}),
         Perm(5, {{1, 2, 4}, {3, 5}}),
         Perm(5, {{1, 2, 4}}),
         Perm(5, {{1, 2, 4, 3}}),
         Perm(5, {{1, 2, 4, 3, 5}}),
         Perm(5, {{1, 2, 4, 5}}),
         Perm(5, {{1, 2, 4, 5, 3}}),
         Perm(5, {{1, 2, 5}, {3, 4}}),
         Perm(5, {{1, 2, 5}}),
         Perm(5, {{1, 2, 5, 3}}),
         Perm(5, {{1, 2, 5, 3, 4}}),
         Perm(5, {{1, 2, 5, 4}}),
         Perm(5, {{1, 2, 5, 4, 3}}),
         Perm(5, {{1, 3}, {2, 4}}),
         Perm(5, {{1, 3}, {2, 4, 5}}),
         Perm(5, {{1, 3}, {2, 5}}),
         Perm(5, {{1, 3}, {2, 5, 4}}),
         Perm(5, {{1, 3}, {4, 5}}),
         Perm(5, {{1, 3}}),
         Perm(5, {{1, 3, 2}, {4,5}}),
         Perm(5, {{1, 3, 2}}),
         Perm(5, {{1, 3, 2, 4}}),
         Perm(5, {{1, 3, 2, 4, 5}}),
         Perm(5, {{1, 3, 2, 5}}),
         Perm(5, {{1, 3, 2, 5, 4}}),
         Perm(5, {{1, 3, 4}, {2, 5}}),
         Perm(5, {{1, 3, 4}}),
         Perm(5, {{1, 3, 4, 2}}),
         Perm(5, {{1, 3, 4, 2, 5}}),
         Perm(5, {{1, 3, 4, 5}}),
         Perm(5, {{1, 3, 4, 5, 2}}),
         Perm(5, {{1, 3, 5}, {2, 4}}),
         Perm(5, {{1, 3, 5}}),
         Perm(5, {{1, 3, 5, 2}}),
         Perm(5, {{1, 3, 5, 2, 4}}),
         Perm(5, {{1, 3, 5, 4}}),
         Perm(5, {{1, 3, 5, 4, 2}}),
         Perm(5, {{1, 4}, {2, 3}}),
         Perm(5, {{1, 4}, {2, 3, 5}}),
         Perm(5, {{1, 4}, {2, 5}}),
         Perm(5, {{1, 4}, {2, 5, 3}}),
         Perm(5, {{1, 4}, {3, 5}}),
         Perm(5, {{1, 4}}),
         Perm(5, {{1, 4, 2}, {3,5}}),
         Perm(5, {{1, 4, 2}}),
         Perm(5, {{1, 4, 2, 3}}),
         Perm(5, {{1, 4, 2, 3, 5}}),
         Perm(5, {{1, 4, 2, 5}}),
         Perm(5, {{1, 4, 2, 5, 3}}),
         Perm(5, {{1, 4, 3}, {2, 5}}),
         Perm(5, {{1, 4, 3}}),
         Perm(5, {{1, 4, 3, 2}}),
         Perm(5, {{1, 4, 3, 2, 5}}),
         Perm(5, {{1, 4, 3, 5}}),
         Perm(5, {{1, 4, 3, 5, 2}}),
         Perm(5, {{1, 4, 5}, {2, 3}}),
         Perm(5, {{1, 4, 5}}),
         Perm(5, {{1, 4, 5, 2}}),
         Perm(5, {{1, 4, 5, 2, 3}}),
         Perm(5, {{1, 4, 5, 3}}),
         Perm(5, {{1, 4, 5, 3, 2}}),
         Perm(5, {{1, 5}, {2, 3}}),
         Perm(5, {{1, 5}, {2, 3, 4}}),
         Perm(5, {{1, 5}, {2, 4}}),
         Perm(5, {{1, 5}, {2, 4, 3}}),
         Perm(5, {{1, 5}, {3, 4}}),
         Perm(5, {{1, 5}}),
         Perm(5, {{1, 5, 2}, {3, 4}}),
         Perm(5, {{1, 5, 2}}),
         Perm(5, {{1, 5, 2, 3}}),
         Perm(5, {{1, 5, 2, 3, 4}}),
         Perm(5, {{1, 5, 2, 4}}),
         Perm(5, {{1, 5, 2, 4, 3}}),
         Perm(5, {{1, 5, 3}, {2, 4}}),
         Perm(5, {{1, 5, 3}}),
         Perm(5, {{1, 5, 3, 2}}),
         Perm(5, {{1, 5, 3, 2, 4}}),
         Perm(5, {{1, 5, 3, 4}}),
         Perm(5, {{1, 5, 3, 4, 2}}),
         Perm(5, {{1, 5, 4}, {2, 3}}),
         Perm(5, {{1, 5, 4}}),
         Perm(5, {{1, 5, 4, 2}}),
         Perm(5, {{1, 5, 4, 2, 3}}),
         Perm(5, {{1, 5, 4, 3}}),
         Perm(5, {{1, 5, 4, 3, 2}}),
         Perm(5, {{2, 3}, {4, 5}}),
         Perm(5, {{2, 3}}),
         Perm(5, {{2, 3, 4}}),
         Perm(5, {{2, 3, 4, 5}}),
         Perm(5, {{2, 3, 5}}),
         Perm(5, {{2, 3, 5, 4}}),
         Perm(5, {{2, 4}, {3, 5}}),
         Perm(5, {{2, 4}}),
         Perm(5, {{2, 4, 3}}),
         Perm(5, {{2, 4, 3, 5}}),
         Perm(5, {{2, 4, 5}}),
         Perm(5, {{2, 4, 5, 3}}),
         Perm(5, {{2, 5}, {3, 4}}),
         Perm(5, {{2, 5}}),
         Perm(5, {{2, 5, 3}}),
         Perm(5, {{2, 5, 3, 4}}),
         Perm(5, {{2, 5, 4}}),
         Perm(5, {{2, 5, 4, 3}}),
         Perm(5, {{3, 4}}),
         Perm(5, {{3, 4, 5}}),
         Perm(5, {{3, 5}}),
         Perm(5, {{3, 5,4}}),
         Perm(5, {{4, 5}})})
    },
    {C1,
      PermutationGroupDescription(1,
        {},
        {Perm(1)})
    },
    {C2,
      PermutationGroupDescription(2,
        {Perm(2, {{1, 2}})},
        {Perm(2),
         Perm(2, {{1, 2}})})
    },
    {C3,
      PermutationGroupDescription(3,
        {Perm(3, {{1, 2, 3}})},
        {Perm(3),
         Perm(3, {{1, 2, 3}}),
         Perm(3, {{1, 3, 2}})})
    },
    {C4,
      PermutationGroupDescription(4,
        {Perm(4, {{1, 2, 3, 4}})},
        {Perm(4),
         Perm(4, {{1, 2, 3, 4}}),
         Perm(4, {{1, 3}, {2, 4}}),
         Perm(4, {{1, 4, 3, 2}})})
    },
    {C5,
      PermutationGroupDescription(5,
        {Perm(5, {{1, 2, 3, 4, 5}})},
        {Perm(5),
         Perm(5, {{1, 2, 3, 4, 5}}),
         Perm(5, {{1, 3, 5, 2, 4}}),
         Perm(5, {{1, 4, 2, 5, 3}}),
         Perm(5, {{1, 5, 4, 3, 2}})})
    },
    {A1,
      PermutationGroupDescription(1,
        {},
        {Perm(1)})
    },
    {A2,
      PermutationGroupDescription(2,
        {},
        {Perm(2)})
    },
    {A3,
      PermutationGroupDescription(3,
        {Perm(3, {{1, 2, 3}})},
        {Perm(3),
         Perm(3, {{1, 2, 3}}),
         Perm(3, {{1, 3, 2}})})
    },
    {A4,
      PermutationGroupDescription(4,
        {Perm(4, {{1, 2, 3}}),
         Perm(4, {{2, 3, 4}})},
        {Perm(4),
         Perm(4, {{1, 2, 3}}),
         Perm(4, {{1, 2, 4}}),
         Perm(4, {{1, 2}, {3, 4}}),
         Perm(4, {{1, 3, 2}}),
         Perm(4, {{1, 3, 4}}),
         Perm(4, {{1, 3}, {2, 4}}),
         Perm(4, {{1, 4, 2}}),
         Perm(4, {{1, 4, 3}}),
         Perm(4, {{1, 4}, {2, 3}}),
         Perm(4, {{2, 3, 4}}),
         Perm(4, {{2, 4, 3}})})
    },
    {A5,
      PermutationGroupDescription(5,
        {Perm(5, {{1, 2, 3, 4, 5}}),
         Perm(5, {{3, 4, 5}})},
        {Perm(5),
         Perm(5, {{1, 2, 3, 4, 5}}),
         Perm(5, {{1, 2, 3, 5, 4}}),
         Perm(5, {{1, 2, 3}}),
         Perm(5, {{1, 2, 4, 3, 5}}),
         Perm(5, {{1, 2, 4, 5, 3}}),
         Perm(5, {{1, 2, 4}}),
         Perm(5, {{1, 2, 5, 3, 4}}),
         Perm(5, {{1, 2, 5, 4, 3}}),
         Perm(5, {{1, 2, 5}}),
         Perm(5, {{1, 2}, {3, 4}}),
         Perm(5, {{1, 2}, {3, 5}}),
         Perm(5, {{1, 2}, {4 ,5}}),
         Perm(5, {{1, 3, 2, 4, 5}}),
         Perm(5, {{1, 3, 2, 5, 4}}),
         Perm(5, {{1, 3, 2}}),
         Perm(5, {{1, 3, 4, 2, 5}}),
         Perm(5, {{1, 3, 4, 5, 2}}),
         Perm(5, {{1, 3, 4}}),
         Perm(5, {{1, 3, 5, 2, 4}}),
         Perm(5, {{1, 3, 5, 4, 2}}),
         Perm(5, {{1, 3, 5}}),
         Perm(5, {{1, 3}, {2, 4}}),
         Perm(5, {{1, 3}, {2, 5}}),
         Perm(5, {{1, 3}, {4, 5}}),
         Perm(5, {{1, 4, 2, 3, 5}}),
         Perm(5, {{1, 4, 2, 5, 3}}),
         Perm(5, {{1, 4, 2}}),
         Perm(5, {{1, 4, 3, 2, 5}}),
         Perm(5, {{1, 4, 3, 5, 2}}),
         Perm(5, {{1, 4, 3}}),
         Perm(5, {{1, 4, 5, 2, 3}}),
         Perm(5, {{1, 4, 5, 3, 2}}),
         Perm(5, {{1, 4, 5}}),
         Perm(5, {{1, 4}, {2, 3}}),
         Perm(5, {{1, 4}, {2, 5}}),
         Perm(5, {{1, 4}, {3, 5}}),
         Perm(5, {{1, 5, 2, 3, 4}}),
         Perm(5, {{1, 5, 2, 4, 3}}),
         Perm(5, {{1, 5, 2}}),
         Perm(5, {{1, 5, 3, 2, 4}}),
         Perm(5, {{1, 5, 3, 4, 2}}),
         Perm(5, {{1, 5, 3}}),
         Perm(5, {{1, 5, 4, 2, 3}}),
         Perm(5, {{1, 5, 4, 3, 2}}),
         Perm(5, {{1, 5, 4}}),
         Perm(5, {{1, 5}, {2, 3}}),
         Perm(5, {{1, 5}, {2, 4}}),
         Perm(5, {{1, 5}, {3, 4}}),
         Perm(5, {{2, 3, 4}}),
         Perm(5, {{2, 3, 5}}),
         Perm(5, {{2, 3}, {4, 5}}),
         Perm(5, {{2, 4}, {3, 5}}),
         Perm(5, {{2, 4, 3}}),
         Perm(5, {{2, 4, 5}}),
         Perm(5, {{2, 5, 3}}),
         Perm(5, {{2, 5, 4}}),
         Perm(5, {{2, 5}, {3, 4}}),
         Perm(5, {{3, 4, 5}}),
         Perm(5, {{3, 5, 4}})})
    },
    {D2,
      PermutationGroupDescription(2,
       {Perm(2, {{1, 2}})},
       {Perm(2),
        Perm(2, {{1, 2}})})
    },
    {D4,
      PermutationGroupDescription(4,
       {Perm(4, {{1, 2}}),
        Perm(4, {{3, 4}})},
       {Perm(4),
        Perm(4, {{1, 2}}),
        Perm(4, {{3, 4}}),
        Perm(4, {{1, 2}, {3, 4}})})
    },
    {D6,
      PermutationGroupDescription(3,
        {Perm(3, {{1, 2, 3}}),
         Perm(3, {{2, 3}})},
        {Perm(3),
         Perm(3, {{1, 2, 3}}),
         Perm(3, {{1, 2}}),
         Perm(3, {{1, 3, 2}}),
         Perm(3, {{1, 3}}),
         Perm(3, {{2, 3}})})
    },
    {D8,
      PermutationGroupDescription(4,
        {Perm(4, {{2, 4}}),
         Perm(4, {{1, 2}, {3, 4}})},
        {Perm(4),
         Perm(4, {{1, 2, 3, 4}}),
         Perm(4, {{1, 2}, {3, 4}}),
         Perm(4, {{1, 3}, {2, 4}}),
         Perm(4, {{1, 3}}),
         Perm(4, {{1, 4, 3, 2}}),
         Perm(4, {{1, 4}, {2, 3}}),
         Perm(4, {{2, 4}})})
    },
    {D10,
      PermutationGroupDescription(5,
        {Perm(5, {{1, 2, 3, 4, 5}}),
         Perm(5, {{2, 5}, {3, 4}})},
        {Perm(5),
         Perm(5, {{1, 5, 4, 3, 2}}),
         Perm(5, {{1, 4, 2, 5, 3}}),
         Perm(5, {{1, 3, 5, 2, 4}}),
         Perm(5, {{1, 2, 3, 4, 5}}),
         Perm(5, {{2, 5}, {3, 4}}),
         Perm(5, {{1, 5}, {2, 4}}),
         Perm(5, {{1, 4}, {2, 3}}),
         Perm(5, {{1, 3}, {4, 5}}),
         Perm(5, {{1, 2}, {3, 5}})})
    },
    {D12,
      PermutationGroupDescription(6,
        {Perm(6, {{1, 2, 3, 4, 5, 6}}),
         Perm(6, {{2, 6}, {3, 5}})},
        {Perm(6),
         Perm(6, {{1, 5, 3}, {2, 6, 4}}),
         Perm(6, {{1, 3, 5}, {2, 4, 6}}),
         Perm(6, {{1, 6, 5, 4, 3, 2}}),
         Perm(6, {{1, 4}, {2, 5}, {3, 6}}),
         Perm(6, {{1, 2, 3, 4, 5, 6}}),
         Perm(6, {{2, 6}, {3, 5}}),
         Perm(6, {{1, 5}, {2, 4}}),
         Perm(6, {{1, 3}, {4, 6}}),
         Perm(6, {{1, 6}, {2, 5}, {3, 4}}),
         Perm(6, {{1, 4}, {2, 3}, {5, 6}}),
         Perm(6, {{1, 2}, {3, 6}, {4, 5}})})
    },
  };

  auto verified_group(verified_groups[group]);

  PermGroup ret(verified_group.degree, verified_group.generators);

#ifndef NDEBUG
  if (!verified_group.verified) {
    assert(perm_group_equal(verified_group.elements, ret)
           && "verified group has correct elements");
    verified_group.verified = true;
  }
#endif

  return ret;
}
