#include <algorithm>
#include <functional>
#include <numeric>
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
#include "util.h"

namespace mpsym
{

TaskMapping ArchGraphSystem::repr_(TaskMapping const &mapping,
                                   TaskOrbits *orbits,
                                   ReprOptions const *options_)
{
  auto options(ReprOptions::fill_defaults(options_));

  auto automs(automorphisms());
  auto gens(automs.generators());

  TaskMapping representative;

  if (options.optimize_symmetric && !_automorphisms_is_shifted_symmetric_valid) {
    _automorphisms_is_shifted_symmetric = automs.is_shifted_symmetric();

    if (_automorphisms_is_shifted_symmetric) {
      _automorphisms_smp = gens.smallest_moved_point();
      _automorphisms_lmp = gens.largest_moved_point();
    }

    _automorphisms_is_shifted_symmetric_valid = true;
  }

  if (options.optimize_symmetric && _automorphisms_is_shifted_symmetric) {
    unsigned task_min = _automorphisms_smp + options.offset;
    unsigned task_max = _automorphisms_lmp + options.offset;

    representative = min_elem_symmetric(mapping, task_min, task_max, &options);

  } else {
    unsigned task_min = 1u + options.offset;
    unsigned task_max = automs.degree() + options.offset;

    representative =
      options.method == ReprMethod::ITERATE ?
        min_elem_iterate(mapping, orbits, &options) :
      options.method == ReprMethod::ORBITS ?
        min_elem_orbits(mapping, orbits, &options) :
      options.method == ReprMethod::LOCAL_SEARCH ?
        options.variant == ReprVariant::LOCAL_SEARCH_SA_LINEAR ?
          min_elem_local_search_sa(mapping, task_min, task_max, &options) :
          min_elem_local_search(mapping, &options) :
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
  auto automs(automorphisms());

  TaskMapping representative(tasks);

  for (auto it = automs.begin(); it != automs.end(); ++it) {
    auto const &factors(it.factors());

    if (tasks.less_than(representative, factors, options->offset))
      representative = tasks.permuted(factors, options->offset);

    if (is_repr(representative, orbits, options)) {
      return representative;
    }
  }

  return representative;
}

TaskMapping ArchGraphSystem::min_elem_orbits(TaskMapping const &tasks,
                                             TaskOrbits *orbits,
                                             ReprOptions const *options)
{
  auto automs(automorphisms());
  auto gens(automs.generators());

  TaskMapping representative(tasks);

  std::unordered_set<TaskMapping> unprocessed, processed;

  unprocessed.insert(tasks);

  while (!unprocessed.empty()) {
    auto it(unprocessed.begin());
    TaskMapping current(*it);
    unprocessed.erase(it);

    processed.insert(current);

    if (current.less_than(representative))
      representative = current;

    for (Perm const &generator : gens) {
      TaskMapping next(current.permuted(generator, options->offset));

      if (is_repr(next, orbits, options))
        return next;
      else if (processed.find(next) == processed.end())
        unprocessed.insert(next);
    }
  }

  return representative;
}

TaskMapping ArchGraphSystem::min_elem_local_search(TaskMapping const &tasks,
                                                   ReprOptions const *options)
{
  auto automs(automorphisms());
  auto gens(local_search_augment_gens(automs, options));

  TaskMapping representative(tasks);

  std::vector<TaskMapping> possible_representatives;
  possible_representatives.reserve(gens.size());

  for (;;) {
    bool stationary = true;

    for (Perm const &generator : gens) {
      if (representative.less_than(representative, generator, options->offset)) {
        if (options->variant == ReprVariant::LOCAL_SEARCH_BFS) {
          possible_representatives.push_back(
            representative.permuted(generator, options->offset));
        } else {
          representative.permute(generator, options->offset);
        }

        stationary = false;
      }
    }

    if (stationary)
      break;

    if (options->variant == ReprVariant::LOCAL_SEARCH_BFS) {
      representative = *std::min_element(possible_representatives.begin(),
                                         possible_representatives.end(),
                                         [](TaskMapping const &lhs,
                                            TaskMapping const &rhs)
                                         { return lhs.less_than(rhs); });

      possible_representatives.clear();
    }
  }

  return representative;
}

PermSet ArchGraphSystem::local_search_augment_gens(
  PermGroup const &automs,
  ReprOptions const *options)
{
  auto gens(automs.generators());

  // append inverse generators
  if (options->local_search_invert_generators) {
    for (auto const &gen : automs.generators())
      gens.insert(~gen);
  }

  // append random generators
  for (unsigned i = 0u; i < options->local_search_append_generators; ++i)
    gens.insert(automs.random_element());

  return gens;
}

TaskMapping ArchGraphSystem::min_elem_local_search_sa(TaskMapping const &tasks,
                                                      unsigned task_min,
                                                      unsigned task_max,
                                                      ReprOptions const *options)
{
  using namespace std::placeholders;

  auto automs(automorphisms());
  auto gens(automs.generators());

  // probability distributions
  static auto re(util::random_engine());

  std::uniform_real_distribution<> d_prob(0.0, 1.0);

  // value function
  auto value(std::bind(local_search_sa_value, _1, task_min, task_max));

  TaskMapping representative(tasks);
  double representative_value = value(representative);

  std::vector<unsigned> gen_indices(gens.size());
  std::iota(gen_indices.begin(), gen_indices.end(), 0u);

  for (unsigned i = 0u; i < options->local_search_sa_iterations; ++i) {
    // schedule T
    double T = local_search_sa_schedule_T(i, options);

    // generate random possible representative
    TaskMapping next;

    auto gen_queue(gen_indices);
    std::shuffle(gen_queue.begin(), gen_queue.end(), re);

    while (!gen_queue.empty()) {
      Perm random_gen(gens[gen_queue.back()]);
      gen_queue.pop_back();

      bool next_valid = false;
      next = representative.permuted(random_gen, options->offset, &next_valid);

      if (next_valid) {
        break;
      }
    }

    // update current representative
    double delta = value(next) - representative_value;

    if (delta > 0.0 || d_prob(re) >= std::exp(delta / T)) {
      representative = next;
      representative_value = value(representative);
    }
  }

  return representative;
}

double ArchGraphSystem::local_search_sa_schedule_T(unsigned i_,
                                                   ReprOptions const *options)
{
  double i = static_cast<double>(i_);
  double i_max = static_cast<double>(options->local_search_sa_iterations);

  double scale = (i_max - i - 1u) / i_max;

  return scale * options->local_search_sa_T_init;
}

double ArchGraphSystem::local_search_sa_value(TaskMapping const &representative,
                                              unsigned task_min,
                                              unsigned task_max)
{
  double ret = 0.0;
  double mult = 1u;

  unsigned num_tasks = 0u;
  for (auto it = representative.rbegin(); it != representative.rend(); ++it) {
    unsigned task = *it;

    if (task < task_min || task > task_max)
      continue;

    ret += mult * (task_max - task);
    mult *= task_max - task_min;

    if (++num_tasks == task_max - task_min + 1u)
      break;
  }

  return std::log(ret - (task_max - task_min)) / num_tasks;
}

TaskMapping ArchGraphSystem::min_elem_symmetric(TaskMapping const &tasks,
                                                unsigned task_min,
                                                unsigned task_max,
                                                ReprOptions const *)
{
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

  return representative;
}

} // namespace mpsym
