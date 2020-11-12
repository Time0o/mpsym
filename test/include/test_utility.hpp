#ifndef _GUARD_TEST_UTILITY_H
#define _GUARD_TEST_UTILITY_H

#include <string>
#include <vector>

#include "gmock/gmock.h"

#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"

enum VerifiedGroup { A1, A2, A3, A4, A5 };

testing::AssertionResult perm_equal(std::vector<unsigned> const &expected,
  mpsym::internal::Perm const &perm);

testing::AssertionResult perm_group_equal(
  mpsym::internal::PermGroup const &expected,
  mpsym::internal::PermGroup const &actual);

testing::AssertionResult perm_group_equal(
  mpsym::internal::PermSet expected_elements,
  mpsym::internal::PermGroup const &actual);

mpsym::internal::PermGroup verified_perm_group(VerifiedGroup group);

#endif // _GUARD_TEST_UTILITY_H
