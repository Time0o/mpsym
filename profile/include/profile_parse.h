#ifndef _GUARD_PROFILE_PARSE_H
#define _GUARD_PROFILE_PARSE_H

#include <memory>
#include <string>
#include <vector>

#include "arch_graph_system.h"
#include "perm.h"
#include "perm_set.h"
#include "permlib.h"
#include "task_allocation.h"

struct GenericGroup
{
  unsigned degree;
  unsigned long long order;
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
    unsigned min_pe;
    unsigned max_pe;
    std::string task_allocations;
  };
}

namespace cgtl
{
  struct TaskAllocationVector
  {
    unsigned min_pe;
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

namespace profile
{

GenericGroup parse_group(std::string const &group_str);

gap::PermSet parse_generators_gap(unsigned degree, std::string const &gen_str);

cgtl::PermSet parse_generators_mpsym(unsigned degree, std::string const &gen_str);

permlib::PermSet parse_generators_permlib(unsigned degree, std::string const &gen_str);

gap::TaskAllocationVector parse_task_allocations_gap(
  std::string const &task_allocations_str);

cgtl::TaskAllocationVector parse_task_allocations_mpsym(
  std::string const &task_allocations_str);

cgtl::TaskAllocationVector parse_task_allocations_gap_to_mpsym(
  std::string const &gap_output_str);

} // namespace profile

#endif // _GUARD_PROFILE_PARSE_H
