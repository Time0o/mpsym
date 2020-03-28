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

namespace gap
{
  using PermGroup = std::string;

  struct PermSet
  {
    unsigned degree;
    std::string permutations;
  };

  using TaskAllocationVector = std::string;
}

namespace mpsym
{
  using TaskAllocationVector = std::vector<mpsym::TaskAllocation>;
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

struct GenericGroup
{
  std::shared_ptr<mpsym::ArchGraphSystem> to_arch_graph_system() const;

  unsigned degree;
  unsigned long long order;
  std::string generators;
};

GenericGroup parse_group(std::string const &group_str);

gap::PermSet parse_generators_gap(unsigned degree, std::string const &gen_str);

mpsym::PermSet parse_generators_mpsym(unsigned degree, std::string const &gen_str);

permlib::PermSet parse_generators_permlib(unsigned degree, std::string const &gen_str);

gap::TaskAllocationVector parse_task_allocations_gap(
  std::string const &task_allocations_str);

mpsym::TaskAllocationVector parse_task_allocations_mpsym(
  std::string const &task_allocations_str);

mpsym::TaskAllocationVector parse_task_allocations_gap_to_mpsym(
  std::vector<std::string> const &gap_output);

} // namespace profile

#endif // _GUARD_PROFILE_PARSE_H
