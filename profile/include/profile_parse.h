#ifndef _GUARD_PROFILE_PARSE_H
#define _GUARD_PROFILE_PARSE_H

#include <string>
#include <vector>

#include "perm.h"
#include "perm_set.h"
#include "permlib.h"
#include "task_mapping.h"

struct GenericGroup
{
  unsigned degree;
  std::string generators;
};

namespace gap
{
  struct PermSet
  {
    unsigned degree;
    std::string permutations;
  };

  struct TaskAllocationVector
  {
    unsigned max_pe;
    std::string task_allocations;
  };
}

namespace cgtl
{
  struct TaskAllocationVector
  {
    unsigned max_pe;
    std::vector<cgtl::TaskAllocation> task_allocations;
  };
}

namespace permlib
{
  struct PermSet
  {
    unsigned degree;
    std::vector<permlib::Permutation::ptr> permutations;
  };
}

GenericGroup parse_group(std::string const &group_str);

gap::PermSet parse_generators_gap(std::string const &gen_str);

cgtl::PermSet parse_generators_mpsym(std::string const &gen_str);

permlib::PermSet parse_generators_permlib(std::string const &gen_str);

gap::TaskAllocationVector parse_task_allocations_gap(
  std::string const &task_allocations_str);

cgtl::TaskAllocationVector parse_task_allocations_mpsym(
  std::string const &task_allocations_str);

#endif // _GUARD_PROFILE_PARSE_H
