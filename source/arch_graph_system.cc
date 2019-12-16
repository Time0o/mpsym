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

PermSet ArchGraphSystem::automorphisms_generators(bool augmented)
{
  static PermSet augmented_generators;
  static bool augmented_generators_valid = false;

  if (!augmented)
    return automorphisms().bsgs().strong_generators();

  if (!augmented_generators_valid) {
    augmented_generators.clear();

    for (Perm const &generator : automorphisms_generators(false)) {
      augmented_generators.insert(generator);
      augmented_generators.insert(~generator);
    }

    augmented_generators_valid = true;
  }

  return augmented_generators;
}

TaskMapping ArchGraphSystem::mapping(TaskMappingRequest const &tmr)
{
  DBG(DEBUG) << "Requested task mapping: " << tmr;

  unsigned min_pe = tmr.offset + 1u;
  unsigned max_pe = min_pe + automorphisms().degree() - 1u;

#ifndef NDEBUG
  if (min_pe != 0u)
    DBG(TRACE) << "Mapping shifted range [" << min_pe << ", " << max_pe << "]";
#endif

  TIMER_CREATE("map approx");
  TIMER_CREATE("map bruteforce");

  TaskAllocation representative =
    tmr.approximate ? min_elem_approx(tmr.allocation, min_pe, max_pe)
                    : min_elem_bruteforce(tmr.allocation, min_pe, max_pe);

  return TaskMapping(tmr.allocation, representative);
}

TaskAllocation ArchGraphSystem::min_elem_bruteforce(TaskAllocation const &tasks,
                                                    unsigned min_pe,
                                                    unsigned max_pe)
{
  DBG(DEBUG) << "Performing brute force mapping";

  TIMER_START("map bruteforce");

  TaskAllocation representative(tasks);

  for (Perm const &element : automorphisms()) {
    if (representative.minimizes(element, min_pe, max_pe))
      representative.permute(element, min_pe, max_pe);
  }

  TIMER_STOP("map bruteforce");

  DBG(DEBUG) << "Found minimal orbit element: " << representative;

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_approx(TaskAllocation const &tasks,
                                                unsigned min_pe,
                                                unsigned max_pe)
{
  DBG(TRACE) << "Performing approximate mapping";

  TaskAllocation representative(tasks);

  TIMER_START("map approx");

  bool stationary = false;
  while (!stationary) {
    stationary = true;

    for (Perm const &generator : automorphisms_generators(true)) {
      if (representative.minimizes(generator, min_pe, max_pe)) {
        representative.permute(generator, min_pe, max_pe);

        stationary = false;
      }
    }
  }

  TIMER_STOP("map approx");

  DBG(DEBUG) << "Found approximate minimal orbit element: " << representative;

  return representative;
}


} // namespace cgtl
