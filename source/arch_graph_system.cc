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

namespace mpsym
{

TaskAllocation ArchGraphSystem::repr(TaskAllocation const &allocation,
                                     ReprOptions const *options_,
                                     TaskOrbits *orbits)
{
  auto options(complete_options(options_));

  DBG(DEBUG) << "Requested task mapping for: " << allocation;

  TaskAllocation representative =
    options.method == ReprMethod::ITERATE ?
      min_elem_iterate(allocation, &options, orbits) :
    options.method == ReprMethod::LOCAL_SEARCH ?
      min_elem_local_search(allocation, &options, orbits) :
    options.method == ReprMethod::ORBITS ?
      min_elem_orbits(allocation, &options, orbits) :
    throw std::logic_error("unreachable");

  if (orbits)
    orbits->insert(representative);

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_iterate(TaskAllocation const &tasks,
                                                 ReprOptions const *options,
                                                 TaskOrbits *orbits)
{
  DBG(DEBUG) << "Performing mapping by iteration";

  TIMER_START("map bruteforce iterate");

  TaskAllocation representative(tasks);

  for (Perm const &element : automorphisms()) {
    if (tasks.less_than(representative, element, options->offset)) {
      representative = tasks.permuted(element, options->offset);

      if (is_repr(representative, options, orbits)) {
        TIMER_STOP("map bruteforce iterate");
        return representative;
      }
    }
  }

  TIMER_STOP("map bruteforce iterate");

  DBG(DEBUG) << "Found minimal orbit element: " << representative;

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_local_search(
  TaskAllocation const &tasks,
  ReprOptions const *options,
  TaskOrbits *)
{
  DBG(TRACE) << "Performing approximate mapping by local search";

  TIMER_START("map approx local search");

  TaskAllocation representative(tasks);

  bool stationary = false;
  while (!stationary) {
    stationary = true;

    for (Perm const &generator : automorphisms().generators()) {
      if (representative.less_than(representative, generator, options->offset)) {
        representative.permute(generator, options->offset);

        stationary = false;
      }
    }
  }

  TIMER_STOP("map approx local search");

  DBG(DEBUG) << "Found approximate minimal orbit element: " << representative;

  return representative;
}

TaskAllocation ArchGraphSystem::min_elem_orbits(TaskAllocation const &tasks,
                                                ReprOptions const *options,
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

    for (Perm const &generator : automorphisms().generators()) {
      TaskAllocation next(current.permuted(generator, options->offset));

      if (is_repr(next, options, orbits)) {
        TIMER_STOP("map bruteforce orbits");
        return next;
      } else if (processed.find(next) == processed.end()) {
        unprocessed.push(next);
      }
    }
  }

  TIMER_STOP("map bruteforce orbits");

  return representative;;
}

} // namespace mpsym
