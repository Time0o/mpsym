#include <algorithm>
#include <queue>
#include <stdexcept>
#include <type_traits>
#include <unordered_set>
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
                                        unsigned offset,
                                        MappingOptions *options,
                                        TaskOrbits *orbits)
{
  options = get_options(options);

  DBG(DEBUG) << "Requested task mapping for: " << allocation;

  TaskAllocation representative =
    options->method == MappingMethod::ITERATE ?
      min_elem_iterate(allocation, offset) :
    options->method == MappingMethod::LOCAL_SEARCH ?
      min_elem_local_search(allocation, offset) :
    options->method == MappingMethod::ORBITS ?
      min_elem_orbits(allocation, offset, orbits) :
    throw std::logic_error("TODO");

  if (orbits)
    orbits->insert(representative);

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_iterate(TaskAllocation const &tasks,
                                                 unsigned offset)
{
  DBG(DEBUG) << "Performing mapping by iteration";

  TIMER_START("map bruteforce iterate");

  TaskAllocation representative(tasks);

  for (Perm const &element : automorphisms()) {
    if (tasks.less_than(representative, element, offset))
      representative = tasks.permuted(element, offset);
  }

  TIMER_STOP("map bruteforce iterate");

  DBG(DEBUG) << "Found minimal orbit element: " << representative;

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_local_search(TaskAllocation const &tasks,
                                                      unsigned offset)
{
  DBG(TRACE) << "Performing approximate mapping by local search";

  TIMER_START("map approx local search");

  TaskAllocation representative(tasks);

  bool stationary = false;
  while (!stationary) {
    stationary = true;

    for (Perm const &generator : automorphisms_generators(true)) {
      if (representative.less_than(representative, generator, offset)) {
        representative.permute(generator, offset);

        stationary = false;
      }
    }
  }

  TIMER_STOP("map approx local search");

  DBG(DEBUG) << "Found approximate minimal orbit element: " << representative;

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_orbits(TaskAllocation const &tasks,
                                                unsigned offset,
                                                TaskOrbits *orbits)
{
  DBG(TRACE) << "Performing mapping by orbit construction";

  TIMER_START("map bruteforce orbits");

  TaskAllocation representative(tasks);

  std::unordered_set<TaskAllocation> processed;
  std::queue<TaskAllocation> unprocessed;

  unprocessed.push(tasks);

  while (!unprocessed.empty()) {
    TaskAllocation current(unprocessed.front());
    unprocessed.pop();

    processed.insert(current);

    if (current.less_than(representative))
      representative = current;

    for (Perm const &generator : automorphisms_generators(false)) {
      TaskAllocation next(current.permuted(generator, offset));

      // TODO
      /*if (std::find(orbits->begin(), orbits->end(), next) != orbits->end()) {
        TIMER_STOP("map bruteforce orbits");
        return next;
      } else*/ if (processed.find(next) == processed.end()) {
        unprocessed.push(next);
      }
    }
  }

  TIMER_STOP("map bruteforce orbits");

  return representative;;
}

ArchGraphSystem::MappingOptions ArchGraphSystem::_default_mapping_options;

} // namespace cgtl
