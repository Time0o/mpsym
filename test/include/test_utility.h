#ifndef _GUARD_TEST_UTILITY_H
#define _GUARD_TEST_UTILITY_H

#include <string>
#include <vector>

#include "gmock/gmock.h"

#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"

enum VerifiedGroup {
  S1, S2, S3, S4, S5,
  C1, C2, C3, C4, C5,
  A1, A2, A3, A4, A5,
  D2, D4, D6, D8, D10,
  D12
};

testing::AssertionResult perm_equal(std::vector<unsigned> const &expected,
  mpsym::Perm const &perm);

testing::AssertionResult perm_group_equal(mpsym::PermGroup const &expected,
                                          mpsym::PermGroup const &actual);

testing::AssertionResult perm_group_equal(mpsym::PermSet expected_elements,
                                          mpsym::PermGroup const &actual);

mpsym::PermGroup verified_perm_group(VerifiedGroup group);

#endif // _GUARD_TEST_UTILITY_H
