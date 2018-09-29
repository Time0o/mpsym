#ifndef _GUARD_TEST_UTILITY_H
#define _GUARD_TEST_UTILITY_H

#include <string>
#include <vector>

#include "gmock/gmock.h"

#include "perm.h"
#include "perm_group.h"

enum VerifiedGroup {
  S1, S2, S3, S4, S5,
  C1, C2, C3, C4, C5,
  A1, A2, A3, A4, A5,
  D1, D2, D3, D4, D5,
  D6, D7, D8, D9, D10
};

testing::AssertionResult perm_equal(std::vector<unsigned> const &expected,
  cgtl::Perm const &perm);

testing::AssertionResult perm_word_equal(std::vector<unsigned> const &expected,
  cgtl::PermWord const &pw);

testing::AssertionResult perm_group_equal(cgtl::PermGroup const &expected,
                                          cgtl::PermGroup const &actual);

testing::AssertionResult perm_group_equal(
  std::vector<std::vector<std::vector<unsigned>>> const &expected,
  cgtl::PermGroup const &actual);

cgtl::PermGroup verified_perm_group(VerifiedGroup group);

std::string resource_path(std::string const &resource);

#endif // _GUARD_TEST_UTILITY_H
