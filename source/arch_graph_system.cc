#include <type_traits>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "arch_graph_system.h"
#include "dbg.h"
#include "perm.h"
#include "perm_set.h"
#include "task_mapping.h"
#include "timer.h"

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

  Timer_create("map approx", Timer::MILLISECONDS);
  Timer_create("map bruteforce", Timer::MILLISECONDS);

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

  Timer_start("map bruteforce");

  TaskAllocation representative(tasks);

  for (Perm const &element : automorphisms()) {
    if (representative.minimizes(element, min_pe, max_pe))
      representative.permute(element, min_pe, max_pe);
  }

  Timer_stop("map bruteforce");

  Dbg(Dbg::DBG) << "Found minimal orbit element: " << representative;

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_approx(TaskAllocation const &tasks,
                                                unsigned min_pe,
                                                unsigned max_pe)
{
  Dbg(Dbg::TRACE) << "Performing approximate mapping";

  TaskAllocation representative(tasks);

  Timer_start("map approx");

  bool stationary = false;
  while (!stationary) {
    stationary = true;

    for (Perm const &generator : automorphisms_generators()) {
      if (representative.minimizes(generator, min_pe, max_pe)) {
        representative.permute(generator, min_pe, max_pe);

        stationary = false;
      }
    }
  }

  Timer_stop("map approx");

  Dbg(Dbg::DBG) << "Found approximate minimal orbit element: " << representative;

  return representative;
}


} // namespace cgtl
