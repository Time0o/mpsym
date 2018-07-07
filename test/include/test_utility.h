#ifndef _GUARD_TEST_UTILITY_H
#define _GUARD_TEST_UTILITY_H

#include <vector>

#include "gmock/gmock.h"
#include "perm.h"
#include "perm_group.h"

::testing::AssertionResult perm_equal(std::vector<unsigned> const &expected,
  cgtl::Perm const &perm);

::testing::AssertionResult perm_word_equal(std::vector<unsigned> const &expected,
  cgtl::PermWord const &pw);

::testing::AssertionResult perm_group_equal(
  std::vector<std::vector<std::vector<unsigned>>> const &expected,
  cgtl::PermGroup const &pg);

#endif // _GUARD_TEST_UTILITY_H
