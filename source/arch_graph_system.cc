#include <type_traits>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "arch_graph_system.h"
#include "dbg.h"
#include "perm.h"
#include "perm_set.h"
#include "task_mapping.h"

namespace
{

using cgtl::Perm;
using cgtl::TaskAllocation;

template<typename AUT>
std::pair<TaskAllocation, bool> find_representative(
  TaskAllocation const &allocation,
  AUT const &automorphisms,
  unsigned min_pe,
  unsigned max_pe,
  bool approximate)
{
  TaskAllocation representative(allocation);
  bool stationary = true;

  for (Perm const &perm : automorphisms) {
    bool new_minimum = false;

    for (auto i = 0u; i < allocation.size(); ++i) {
      unsigned task_base = allocation[i];
      if (task_base < min_pe || task_base > max_pe)
        continue;

      unsigned task = task_base - (min_pe - 1u);

      unsigned permuted = perm[task] + (min_pe - 1u);

      if (permuted < representative[i]) {
        new_minimum = true;
        break;
      }

      if (permuted > representative[i])
        break;
    }

    if (new_minimum) {
      for (auto i = 0u; i < allocation.size(); ++i) {
        unsigned task_base = allocation[i];
        if (task_base < min_pe || task_base > max_pe)
          continue;

        unsigned task = task_base - (min_pe - 1u);

        representative[i] = perm[task] + (min_pe - 1u);
      }

      stationary = false;

      if (approximate)
        break;
    }
  }

  return {representative, stationary};
}

} // namespace

namespace cgtl
{

TaskMapping ArchGraphSystem::mapping(TaskMappingRequest const &tmr)
{
  Dbg(Dbg::DBG) << "Requested task mapping: " << tmr;

  unsigned min_pe = tmr.offset + 1u;
  unsigned max_pe = min_pe + automorphisms().degree() - 1u;

#ifndef NDEBUG
  if (min_pe != 0u) {
    Dbg(Dbg::TRACE) << "Mapping shifted range ["
                    << min_pe << ", " << max_pe << "]";
  }
#endif

  TaskAllocation representative =
    tmr.approximate ? min_elem_approx(tmr.allocation, min_pe, max_pe)
                    : min_elem_bruteforce(tmr.allocation, min_pe, max_pe);

  return TaskMapping(tmr.allocation, representative);
}

TaskAllocation ArchGraphSystem::min_elem_bruteforce(TaskAllocation const &tasks,
                                                    unsigned min_pe,
                                                    unsigned max_pe)
{

  Dbg(Dbg::DBG) << "Performing brute force mapping";

  auto repr(find_representative(tasks, automorphisms(), min_pe, max_pe, false));
  TaskAllocation representative = repr.first;

  Dbg(Dbg::DBG) << "Found minimal orbit element: " << representative;

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_approx(TaskAllocation const &tasks,
                                                unsigned min_pe,
                                                unsigned max_pe)
{
  Dbg(Dbg::TRACE) << "Performing approximate mapping";

  TaskAllocation representative(tasks);
  bool stationary;

  for (;;) {
    auto repr(find_representative(
      representative, automorphisms_generators(), min_pe, max_pe, true));

    representative = repr.first;
    stationary = repr.second;

    if (stationary)
      break;
  }

  Dbg(Dbg::DBG) << "Found approximate minimal orbit element: " << representative;

  return representative;
}


} // namespace cgtl
