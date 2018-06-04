#ifndef _GUARD_TEST_UTILITY_H
#define _GUARD_TEST_UTILITY_H

#include <vector>

#include "gmock/gmock.h"
#include "perm.h"

::testing::AssertionResult perm_equal(std::vector<unsigned> const &expected,
  cgtl::Perm const &perm);

#endif // _GUARD_TEST_UTILITY_H
