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
#include "task_mapping.h"
#include "task_orbits.h"
#include "timer.h"

namespace mpsym
{

TaskMapping ArchGraphSystem::repr_(TaskMapping const &mapping,
                                   TaskOrbits *orbits,
                                   ReprOptions const *options_)
{
  auto options(ReprOptions::fill_defaults(options_));

  TaskMapping representative;

  if (options.optimize_symmetric && automorphisms().is_shifted_symmetric()) {
    auto generators(automorphisms().generators());

    unsigned task_min = generators.smallest_moved_point() + options.offset;
    unsigned task_max = generators.largest_moved_point() + options.offset;

    representative = min_elem_symmetric(mapping, task_min, task_max, &options);

  } else {
    representative = options.method == ReprMethod::ITERATE ?
                       min_elem_iterate(mapping, orbits, &options) :
                     options.method == ReprMethod::LOCAL_SEARCH ?
                       min_elem_local_search(mapping, &options) :
                     options.method == ReprMethod::ORBITS ?
                       min_elem_orbits(mapping, orbits, &options) :
                     throw std::logic_error("unreachable");
  }

  if (orbits)
    orbits->insert(representative);

  return representative;
}

TaskMapping ArchGraphSystem::min_elem_iterate(TaskMapping const &tasks,
                                              TaskOrbits *orbits,
                                              ReprOptions const *options)
{
  TIMER_START("map bruteforce iterate");

  TaskMapping representative(tasks);

  for (Perm const &element : automorphisms()) {
    if (tasks.less_than(representative, element, options->offset)) {
      representative = tasks.permuted(element, options->offset);

      if (is_repr(representative, orbits, options)) {
        TIMER_STOP("map bruteforce iterate");
        return representative;
      }
    }
  }

  TIMER_STOP("map bruteforce iterate");

  return representative;
}

TaskMapping ArchGraphSystem::min_elem_orbits(TaskMapping const &tasks,
                                             TaskOrbits *orbits,
                                             ReprOptions const *options)
{
  TIMER_START("map bruteforce orbits");

  TaskMapping representative(tasks);

  std::unordered_set<TaskMapping> processed;
  std::queue<TaskMapping> unprocessed;

  unprocessed.push(tasks);

  while (!unprocessed.empty()) {
    TaskMapping current(unprocessed.front());
    unprocessed.pop();

    processed.insert(current);

    if (current.less_than(representative))
      representative = current;

    for (Perm const &generator : automorphisms().generators()) {
      TaskMapping next(current.permuted(generator, options->offset));

      if (is_repr(next, orbits, options)) {
        TIMER_STOP("map bruteforce orbits");
        return next;
      } else if (processed.find(next) == processed.end()) {
        unprocessed.push(next);
      }
    }
  }

  TIMER_STOP("map bruteforce orbits");

  return representative;
}

TaskMapping ArchGraphSystem::min_elem_local_search(TaskMapping const &tasks,
                                                   ReprOptions const *options)
{
  TIMER_START("map approx local search");

  TaskMapping representative(tasks);

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

  return representative;
}

TaskMapping ArchGraphSystem::min_elem_symmetric(TaskMapping const &tasks,
                                                unsigned task_min,
                                                unsigned task_max,
                                                ReprOptions const *)
{
  TIMER_START("map symmetric");

  TaskMapping representative(tasks);

  std::vector<unsigned> perm(task_max - task_min + 1u);
  unsigned perm_next = task_min;

  for (auto i = 0u; i < tasks.size(); ++i) {
    unsigned task = tasks[i];
    if (task < task_min || task > task_max)
      continue;

    auto j = task - task_min;

    unsigned to = perm[j];

    if (to == 0u) {
      to = perm[j] = perm_next;
      ++perm_next;
    }

    representative[i] = to;
  }

  TIMER_STOP("map symmetric");

  return representative;
}

} // namespace mpsym
