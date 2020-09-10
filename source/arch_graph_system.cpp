#include <algorithm>
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

unsigned ArchGraphSystem::num_automorphism_orbits(
  unsigned num_tasks,
  bool unique_tasks,
  AutomorphismOptions const *options)
{
  auto automs(automorphisms(options));

  BSGS::order_type ret = 0;

  for (Perm const &perm : automs) {
    unsigned pes_fixed = 0u;
    for (unsigned pe = 1u; pe <= automs.degree(); ++pe) {
      if (perm[pe] == pe)
        ++pes_fixed;
    }

    if (unique_tasks) {
      if (pes_fixed < num_tasks)
        continue;

      BSGS::order_type task_mappings_moved = 1;
      for (unsigned n = pes_fixed; n > pes_fixed - num_tasks; --n)
        task_mappings_moved *= n;

      ret += task_mappings_moved;

    } else {
      ret += pow(static_cast<BSGS::order_type>(pes_fixed), num_tasks);
    }
  }

  ret /= automs.order();

  if (ret > std::numeric_limits<unsigned>::max())
    throw std::runtime_error("orbit size limit reached");

  return static_cast<unsigned>(ret);
}

std::vector<unsigned> ArchGraphSystem::automorphism_orbit_sizes(
  unsigned num_tasks,
  bool unique_tasks,
  AutomorphismOptions const *options)
{
  auto automs(automorphisms(options));

  // number of orbits
  auto num_orbits(num_automorphism_orbits(num_tasks, unique_tasks, options));

  // number of processing element indices stabilized by all automorphisms
  std::vector<int> stabilized(automs.degree(), 1);
  for (Perm const &perm : automs) {
    for (unsigned p = 1u; p <= automs.degree(); ++p) {
      if (stabilized[p - 1u] && perm[p] != p)
        stabilized[p - 1u] = 0;
    }
  }

  unsigned num_stabilized =
    std::accumulate(stabilized.begin(), stabilized.end(), 0u);

  // number of singular orbits
  BSGS::order_type num_singular_orbits = 0;

  if (unique_tasks) {
    if (num_stabilized > num_tasks) {
      num_singular_orbits = 1;
      for (unsigned n = num_stabilized; n > num_stabilized - num_tasks; --n)
        num_singular_orbits *= n;
    }
  } else {
    num_singular_orbits = pow(static_cast<BSGS::order_type>(num_stabilized), num_tasks);
  }

  // find all non singular orbits
  std::vector<std::unordered_set<TaskMapping>> non_singular_orbits;

  try {
    while (non_singular_orbits.size() < num_orbits - num_singular_orbits) {
      TaskMapping next_representative;

      for (;;) {
        next_representative = random_task_mapping(num_tasks, unique_tasks);

        bool new_orbit = true;
        for (auto const &orbit : non_singular_orbits) {
          if (orbit.find(next_representative) != orbit.end()) {
            new_orbit = false;
            break;
          }
        }

        if (new_orbit)
          break;
      }

      auto orbit(task_mapping_orbit(next_representative));

      if (orbit.size() > 1u)
        non_singular_orbits.emplace_back(std::move(orbit));
    }
  } catch (std::bad_alloc const &) {
    non_singular_orbits.clear();
    non_singular_orbits.shrink_to_fit();

    return {};
  }

  std::vector<unsigned> ret(num_orbits, 1u);
  for (auto i = 0u; i < non_singular_orbits.size(); ++i)
    ret[i] = non_singular_orbits[i].size();

  std::sort(ret.begin(), ret.end());

  return ret;
}

std::vector<TaskMapping> ArchGraphSystem::orbit(TaskMapping const &mapping)
{
  auto automs(automorphisms());
  auto gens(automs.generators());

  std::unordered_set<TaskMapping> unprocessed, processed;

  unprocessed.insert(mapping);

  while (!unprocessed.empty()) {
    auto it(unprocessed.begin());
    TaskMapping current(*it);
    unprocessed.erase(it);

    processed.insert(current);

    for (Perm const &gen : gens) {
      TaskMapping next(current.permuted(gen));

      if (processed.find(next) == processed.end())
        unprocessed.insert(next);
    }
  }

  return std::vector<TaskMapping>(processed.begin(), processed.end());
}

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
      options.method == ReprOptions::Method::ITERATE ?
        min_elem_iterate(mapping, orbits, &options) :
      options.method == ReprOptions::Method::ORBITS ?
        min_elem_orbits(mapping, orbits, &options) :
      options.method == ReprOptions::Method::LOCAL_SEARCH ?
        options.variant == ReprOptions::Variant::LOCAL_SEARCH_SA_LINEAR ?
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

    for (Perm const &gen : gens) {
      TaskMapping next(current.permuted(gen, options->offset));

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

    for (Perm const &gen : gens) {
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

TaskMapping ArchGraphSystem::random_task_mapping(unsigned num_tasks,
                                                 bool unique_tasks)
{
  static auto re(util::random_engine());

  auto automs(automorphisms());

  std::uniform_int_distribution<unsigned> d(1u, automs.degree());

  std::vector<unsigned> tasks(num_tasks);
  std::vector<int> included(automs.degree(), 0);

  for (unsigned &task : tasks) {
    do {
      task = d(re);
    } while (unique_tasks && included[task - 1u]);

    included[task - 1u] = 1;
  }

  return TaskMapping(tasks);
}

std::unordered_set<TaskMapping> ArchGraphSystem::task_mapping_orbit(
  TaskMapping const &tasks)
{
  auto automs(automorphisms());
  auto gens(automs.generators());

  std::unordered_set<TaskMapping> unprocessed, orbit;

  unprocessed.insert(tasks);

  while (!unprocessed.empty()) {
    auto it(unprocessed.begin());
    TaskMapping current(*it);
    unprocessed.erase(it);

    if (orbit.size() == std::numeric_limits<unsigned>::max())
      throw std::runtime_error("orbit size limit reached");

    orbit.insert(current);

    for (Perm const &generator : gens) {
      TaskMapping next(current.permuted(generator));

      if (orbit.find(next) == orbit.end())
        unprocessed.insert(next);
    }
  }

  return orbit;
}

} // namespace mpsym
