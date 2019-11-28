#ifndef _GUARD_PROFILE_PARSE_H
#define _GUARD_PROFILE_PARSE_H

#include <string>
#include <tuple>
#include <vector>

#include "perm.h"
#include "perm_set.h"
#include "permlib.h"
#include "task_mapping.h"

namespace permlib
{
  struct PermSet
  {
    unsigned degree;
    std::vector<permlib::Permutation::ptr> permutations;
  };
}

std::tuple<unsigned, unsigned, std::string> parse_group(
  std::string const &group_str);

std::string parse_generators_gap(std::string const &gen_str);

cgtl::PermSet parse_generators_mpsym(std::string const &gen_str);

permlib::PermSet parse_generators_permlib(std::string const &gen_str);

std::string parse_task_allocations_gap(std::string const &task_allocations_str);

std::vector<cgtl::TaskAllocation> parse_task_allocations_mpsym(
  std::string const &task_allocations_str);

#endif // _GUARD_PROFILE_PARSE_H
