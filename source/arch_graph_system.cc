#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "arch_graph_system.h"
#include "dbg.h"
#include "perm.h"
#include "perm_set.h"
#include "task_mapping.h"
#include "task_orbits.h"
#include "timer.h"

namespace cgtl
{

PermSet ArchGraphSystem::automorphisms_generators(bool augmented)
{
  if (!augmented)
    return automorphisms().bsgs().strong_generators();

  if (!_augmented_generators_valid) {
    _augmented_generators.clear();

    for (Perm const &generator : automorphisms_generators(false)) {
      _augmented_generators.insert(generator);

      Perm generator_inverted(~generator);
      if (generator_inverted != generator)
        _augmented_generators.insert(generator_inverted);
    }

    _augmented_generators_valid = true;
  }

  return _augmented_generators;
}

TaskMapping ArchGraphSystem::mapping(TaskMappingRequest const &tmr,
                                     MappingMethod method,
                                     TaskOrbits *orbits)
{
  DBG(DEBUG) << "Requested task mapping: " << tmr;

  unsigned min_pe = tmr.offset + 1u;
  unsigned max_pe = min_pe + automorphisms().degree() - 1u;

#ifndef NDEBUG
  if (min_pe != 0u) {
    DBG(TRACE) << "Mapping shifted range [" << min_pe << ", " << max_pe << "]";
  }
#endif

  TaskAllocation representative =
    method == MappingMethod::BRUTEFORCE ?
      min_elem_bruteforce(tmr.allocation, min_pe, max_pe) :
    method == MappingMethod::APPROXIMATE ?
      min_elem_approx(tmr.allocation, min_pe, max_pe) :
    throw std::logic_error("unreachable");

  if (orbits)
    orbits->insert(representative);

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
    bool new_minimum = tasks.permutes_to_less_than(
      representative, {element, min_pe, max_pe});

    if (new_minimum)
      representative = tasks.permuted({element, min_pe, max_pe});
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
      bool new_minimum = representative.permutes_to_less_than(
        representative, {generator, min_pe, max_pe});

      if (new_minimum) {
        representative.permute({generator, min_pe, max_pe});

        stationary = false;
      }
    }
  }

  TIMER_STOP("map approx");

  DBG(DEBUG) << "Found approximate minimal orbit element: " << representative;

  return representative;
}


} // namespace cgtl
