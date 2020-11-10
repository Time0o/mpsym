#include <algorithm>
#include <atomic>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "arch_graph.hpp"
#include "arch_graph_automorphisms.hpp"
#include "arch_graph_cluster.hpp"
#include "arch_graph_system.hpp"
#include "arch_uniform_super_graph.hpp"
#include "bsgs.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"
#include "task_mapping.hpp"
#include "task_orbits.hpp"
#include "timeout.hpp"
#include "util.hpp"

using boost::multiprecision::pow;

namespace mpsym
{

using namespace internal;

std::shared_ptr<ArchGraphSystem> ArchGraphSystem::expand_automorphisms() const
{
  auto const *ag(dynamic_cast<ArchGraph const *>(this));
  if (ag) {
    auto ag_copy(std::make_shared<ArchGraph>(*ag));
    return std::make_shared<ArchGraphAutomorphisms>(ag_copy->automorphisms());
  }

  auto const *agc(dynamic_cast<ArchGraphCluster const *>(this));
  if (agc) {
    auto agc_copy(std::make_shared<ArchGraphCluster>());
    for (auto ss : agc->subsystems())
      agc_copy->add_subsystem(ss->expand_automorphisms());

    return agc_copy;
  }

  auto const *ausg(dynamic_cast<ArchUniformSuperGraph const *>(this));
  if (ausg) {
    auto ausg_copy(std::make_shared<ArchUniformSuperGraph>(
      ausg->super_graph()->expand_automorphisms(),
      ausg->proto()->expand_automorphisms()));

    return ausg_copy;
  }

  auto const *aga(dynamic_cast<ArchGraphAutomorphisms const *>(this));
  if (aga) {
    auto aga_copy(std::make_shared<ArchGraphAutomorphisms>(*aga));
    return aga_copy;
  }

  throw std::logic_error("unreachable");
}

std::vector<TaskMapping> ArchGraphSystem::orbit(TaskMapping const &mapping)
{
  automorphisms();

  std::unordered_set<TaskMapping> unprocessed, processed;

  unprocessed.insert(mapping);

  while (!unprocessed.empty()) {
    auto it(unprocessed.begin());
    TaskMapping current(*it);
    unprocessed.erase(it);

    processed.insert(current);

    for (Perm const &gen : _automorphism_generators) {
      TaskMapping next(current.permuted(gen));

      if (processed.find(next) == processed.end())
        unprocessed.insert(next);
    }
  }

  return std::vector<TaskMapping>(processed.begin(), processed.end());
}

bool ArchGraphSystem::automorphisms_symmetric(ReprOptions const *options)
{
  TaskMapping representative;

  if (options->optimize_symmetric && !_automorphisms_is_shifted_symmetric_valid) {
    _automorphisms_is_shifted_symmetric = _automorphisms.is_shifted_symmetric();

    if (_automorphisms_is_shifted_symmetric) {
      _automorphisms_smp = _automorphism_generators.smallest_moved_point();
      _automorphisms_lmp = _automorphism_generators.largest_moved_point();
    }

    _automorphisms_is_shifted_symmetric_valid = true;
  }

  return options->optimize_symmetric && _automorphisms_is_shifted_symmetric;
}

TaskMapping ArchGraphSystem::repr_(TaskMapping const &mapping,
                                   ReprOptions const *options_,
                                   TaskOrbits *orbits,
                                   timeout::flag aborted)
{
  automorphisms();

  auto options(ReprOptions::fill_defaults(options_));

  if (_automorphisms.is_trivial())
    return mapping;

  if (automorphisms_symmetric(&options))
    return min_elem_symmetric(mapping, &options);

  return options.method == ReprOptions::Method::ITERATE ?
           min_elem_iterate(mapping, &options, orbits, aborted) :
         options.method == ReprOptions::Method::ORBITS ?
           min_elem_orbits(mapping, &options, orbits, aborted) :
         options.method == ReprOptions::Method::LOCAL_SEARCH ?
           options.variant == ReprOptions::Variant::LOCAL_SEARCH_SA_LINEAR ?
             min_elem_local_search_sa(mapping, &options) :
             min_elem_local_search(mapping, &options) :
         throw std::logic_error("unreachable");
}

TaskMapping ArchGraphSystem::min_elem_iterate(TaskMapping const &tasks,
                                              ReprOptions const *options,
                                              TaskOrbits *orbits,
                                              timeout::flag aborted) const
{
  TaskMapping representative(tasks);

  for (auto it = _automorphisms.begin(); it != _automorphisms.end(); ++it) {
    if (timeout::is_set(aborted))
      throw timeout::AbortedError("min_elem_iterate");

    auto const &factors(it.factors());

    if (tasks.less_than(representative, factors, options->offset))
      representative = tasks.permuted(factors, options->offset);

    if (is_repr(representative, options, orbits)) {
      return representative;
    }
  }

  return representative;
}

TaskMapping ArchGraphSystem::min_elem_orbits(TaskMapping const &tasks,
                                             ReprOptions const *options,
                                             TaskOrbits *orbits,
                                             timeout::flag aborted) const
{
  TaskMapping representative(tasks);

  std::unordered_set<TaskMapping> unprocessed, processed;

  unprocessed.insert(tasks);

  while (!unprocessed.empty()) {
    if (timeout::is_set(aborted))
      throw timeout::AbortedError("min_elem_orbits");

    auto it(unprocessed.begin());
    TaskMapping current(*it);
    unprocessed.erase(it);

    processed.insert(current);

    if (current.less_than(representative))
      representative = current;

    for (Perm const &gen : _automorphism_generators) {
      TaskMapping next(current.permuted(gen, options->offset));

      if (is_repr(next, options, orbits))
        return next;
      else if (processed.find(next) == processed.end())
        unprocessed.insert(next);
    }
  }

  return representative;
}

TaskMapping ArchGraphSystem::min_elem_local_search(
  TaskMapping const &tasks,
  ReprOptions const *options) const
{
  auto generators(local_search_augment_gens(options));

  TaskMapping representative(tasks);

  std::vector<TaskMapping> possible_representatives;
  possible_representatives.reserve(generators.size());

  for (;;) {
    bool stationary = true;

    for (Perm const &gen : generators) {
      if (representative.less_than(representative, gen, options->offset)) {
        if (options->variant == ReprOptions::Variant::LOCAL_SEARCH_BFS) {
          possible_representatives.push_back(
            representative.permuted(gen, options->offset));
        } else {
          representative.permute(gen, options->offset);
        }

        stationary = false;
      }
    }

    if (stationary)
      break;

    if (options->variant == ReprOptions::Variant::LOCAL_SEARCH_BFS) {
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
  ReprOptions const *options) const
{
  auto generators(_automorphism_generators);

  // append inverse generators
  if (options->local_search_invert_generators) {
    for (auto const &gen : _automorphisms.generators())
      generators.insert(~gen);
  }

  // append random generators
  for (unsigned i = 0u; i < options->local_search_append_generators; ++i)
    generators.insert(_automorphisms.random_element());

  return generators;
}

TaskMapping ArchGraphSystem::min_elem_local_search_sa(
  TaskMapping const &tasks,
  ReprOptions const *options) const
{
  using namespace std::placeholders;

  // probability distributions
  static auto re(util::random_engine());

  std::uniform_real_distribution<> d_prob(0.0, 1.0);

  // value function
  unsigned task_min = 1u + options->offset;
  unsigned task_max = _automorphisms.degree() + options->offset;

  auto value(std::bind(local_search_sa_value, _1, task_min, task_max));

  TaskMapping representative(tasks);
  double representative_value = value(representative);

  std::vector<unsigned> gen_indices(_automorphism_generators.size());
  std::iota(gen_indices.begin(), gen_indices.end(), 0u);

  for (unsigned i = 0u; i < options->local_search_sa_iterations; ++i) {
    // schedule T
    double T = local_search_sa_schedule_T(i, options);

    // generate random possible representative
    TaskMapping next;

    auto gen_queue(gen_indices);
    std::shuffle(gen_queue.begin(), gen_queue.end(), re);

    while (!gen_queue.empty()) {
      Perm random_gen(_automorphism_generators[gen_queue.back()]);
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

TaskMapping ArchGraphSystem::min_elem_symmetric(
  TaskMapping const &tasks,
  ReprOptions const *options) const
{
  TaskMapping representative(tasks);

  unsigned task_min = _automorphisms_smp + options->offset;
  unsigned task_max = _automorphisms_lmp + options->offset;

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
