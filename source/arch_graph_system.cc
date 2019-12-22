#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "arch_graph_system.h"
#include "dbg.h"
#include "perm.h"
#include "perm_set.h"
#include "task_allocation.h"
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

TaskAllocation ArchGraphSystem::mapping(TaskAllocation const &allocation,
                                        MappingMethod method,
                                        TaskOrbits *orbits)
{
  DBG(DEBUG) << "Requested task mapping for: " << allocation;

  TaskAllocation representative =
    method == MappingMethod::BRUTEFORCE ? min_elem_bruteforce(allocation) :
    method == MappingMethod::APPROXIMATE ? min_elem_approx(allocation) :
    throw std::logic_error("unreachable");

  if (orbits)
    orbits->insert(representative);

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_bruteforce(TaskAllocation const &tasks)
{
  DBG(DEBUG) << "Performing brute force mapping";

  TIMER_START("map bruteforce");

  TaskAllocation representative(tasks);

  for (Perm const &element : automorphisms()) {
    if (tasks.less_than(representative, element))
      representative = tasks.permuted(element);
  }

  TIMER_STOP("map bruteforce");

  DBG(DEBUG) << "Found minimal orbit element: " << representative;

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_approx(TaskAllocation const &tasks)
{
  DBG(TRACE) << "Performing approximate mapping";

  TIMER_START("map approx");

  TaskAllocation representative(tasks);

  bool stationary = false;
  while (!stationary) {
    stationary = true;

    for (Perm const &generator : automorphisms_generators(true)) {
      if (representative.less_than(representative, generator)) {
        representative.permute(generator);

        stationary = false;
      }
    }
  }

  TIMER_STOP("map approx");

  DBG(DEBUG) << "Found approximate minimal orbit element: " << representative;

  return representative;
}


} // namespace cgtl
