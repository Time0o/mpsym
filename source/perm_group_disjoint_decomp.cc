#include <algorithm>
#include <tuple>
#include <vector>

#include "dbg.h"
#include "orbits.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"
#include "timer.h"

/**
 * @file perm_group_disjoint_decomp.cc
 * @brief Implements `PermGroup::disjoint_decomposition`.
 *
 * @author Timo Nicolai
 */

namespace mpsym
{

namespace internal
{

std::vector<PermGroup> PermGroup::disjoint_decomposition(
  bool complete,
  bool disjoint_orbit_optimization) const
{
  return complete ? disjoint_decomp_complete(disjoint_orbit_optimization)
                  : disjoint_decomp_incomplete();
}

bool PermGroup::disjoint_decomp_orbits_dependent(
  Orbit const &orbit1,
  Orbit const &orbit2) const
{
  std::unordered_set<Perm> restricted_stabilizers, restricted_elements;

  for (Perm const &perm : *this) {
    Perm restricted_perm(perm.restricted(orbit1.begin(), orbit1.end()));

    if (restricted_perm.id())
      continue;

    if (perm.stabilizes(orbit2.begin(), orbit2.end()))
      restricted_stabilizers.insert(restricted_perm);

    restricted_elements.insert(restricted_perm);
  }

  return restricted_stabilizers.size() < restricted_elements.size();
}

void PermGroup::disjoint_decomp_generate_dependency_classes(
  OrbitPartition &orbits) const
{
  unsigned num_dependency_classes = 0u;

  std::vector<int> processed(orbits.num_partitions(), 0);
  unsigned num_processed = 0u;

  for (unsigned i = 0u; i < orbits.num_partitions(); ++i) {
    if (processed[i])
      continue;

    // determine which orbits to merge
    std::unordered_set<unsigned> merge{i};

    for (unsigned j = i + 1u; j < orbits.num_partitions(); ++j) {
      if (processed[j])
        continue;

      if (disjoint_decomp_orbits_dependent(orbits[i], orbits[j])) {
        merge.insert(j);

        processed[j] = 1;
        ++num_processed;
      }
    }

    // merge orbits
    for (unsigned x = 1u; x <= degree(); ++x) {
      if (merge.find(orbits.partition_index(x)) != merge.end())
        orbits.change_partition(x, num_dependency_classes);
    }

    ++num_dependency_classes;

    // check if we're done
    processed[i] = 1;

    if (++num_processed == orbits.num_partitions())
      break;
  }
}

bool PermGroup::disjoint_decomp_restricted_subgroups(
  OrbitPartition const &orbit_split,
  PermGroup const &perm_group,
  std::pair<PermGroup, PermGroup> &restricted_subgroups)
{
  auto split1(orbit_split[0]);
  auto split2(orbit_split[1]);

  PermSet restricted_generators1;
  PermSet restricted_generators2;

  for (Perm const &gen : perm_group.generators()) {
    Perm restricted_generator1(gen.restricted(split1.begin(), split1.end()));
    Perm restricted_generator2(gen.restricted(split2.begin(), split2.end()));

    if (!perm_group.contains_element(restricted_generator1) ||
        !perm_group.contains_element(restricted_generator2)) {

      DBG(TRACE) << "Restricted groups are not a disjoint subgroup decomposition";

      return false;
    }

    restricted_generators1.insert(restricted_generator1);
    restricted_generators2.insert(restricted_generator2);
  }

  restricted_subgroups.first = PermGroup(perm_group.degree(),
                                         restricted_generators1);

  restricted_subgroups.second = PermGroup(perm_group.degree(),
                                          restricted_generators2);

  DBG(TRACE) << "Found disjoint subgroup decomposition:";
  DBG(TRACE) << restricted_subgroups.first;
  DBG(TRACE) << restricted_subgroups.second;

  return true;
}

std::vector<PermGroup> PermGroup::disjoint_decomp_join_results(
  std::vector<PermGroup> const &res1,
  std::vector<PermGroup> const &res2)
{
  auto res(res1);

  res.insert(res.end(), res2.begin(), res2.end());

  return res;
}

std::vector<PermGroup> PermGroup::disjoint_decomp_complete_recursive(
  OrbitPartition const &orbits,
  PermGroup const &perm_group)
{
  // iterate over all possible partitions of the set of all orbits into two sets
  assert(orbits.num_partitions() < 8 * sizeof(unsigned long long));

  for (auto part = 1ULL; !(part & (1ULL << (orbits.num_partitions() - 1u))); ++part) {
    OrbitPartition orbit_split(perm_group.degree());

    for (unsigned x = 1u; x <= perm_group.degree(); ++x) {
      if (orbits.partition_index(x) == -1)
        continue;

      if ((1ULL << orbits.partition_index(x)) & part)
        orbit_split.change_partition(x, 1);
      else
        orbit_split.change_partition(x, 0);
    }

    DBG(TRACE) << "Considering orbit split:";
    DBG(TRACE) << orbit_split;

    // try to find restricted subgroup decomposition
    std::pair<PermGroup, PermGroup> restricted_subgroups;

    if (!disjoint_decomp_restricted_subgroups(
          orbit_split, perm_group, restricted_subgroups))
      continue;

    DBG(TRACE) << "Restricted groups are a disjoint subgroup decomposition";

    // recurse for both orbit partition elements and return combined result
    auto orbits_recurse(orbits.split(orbit_split));

    DBG(TRACE) << "Recursing with orbit partitions:";
    DBG(TRACE) << orbits_recurse[0];
    DBG(TRACE) << orbits_recurse[1];

    return disjoint_decomp_join_results(
      disjoint_decomp_complete_recursive(orbits_recurse[0],
                                         restricted_subgroups.first),
      disjoint_decomp_complete_recursive(orbits_recurse[1],
                                         restricted_subgroups.second));
  }

  DBG(TRACE) << "No further decomposition possible, returning group";

  return {perm_group};
}

std::vector<PermGroup> PermGroup::disjoint_decomp_complete(
  bool disjoint_orbit_optimization) const
{
  DBG(DEBUG) << "Finding (complete) disjoint subgroup decomposition for:";
  DBG(DEBUG) << *this;

  OrbitPartition orbits(degree(), generators());

  DBG(TRACE) << "Orbit decomposition:";
  DBG(TRACE) << orbits;

  if (disjoint_orbit_optimization) {
    DBG(TRACE) << "Using dependent orbit optimization";
    disjoint_decomp_generate_dependency_classes(orbits);

    DBG(TRACE) << "=> Grouped dependency class unions:";
    DBG(TRACE) << orbits;
  }

  auto decomp(disjoint_decomp_complete_recursive(orbits, *this));

  DBG(DEBUG) << "Found disjoint subgroup decomposition:";
  for (PermGroup const &pg : decomp)
    DBG(DEBUG) << pg;

  return decomp;
}

void
PermGroup::MovedSet::init(Perm const &perm)
{
  clear();
  for (unsigned i = 1u; i <= perm.degree(); ++i) {
    if (perm[i] != i)
      push_back(i);
  }
}

bool PermGroup::MovedSet::equivalent(MovedSet const &other) const
{
  MovedSet::size_type i1 = 0u, i2 = 0u;

  while ((*this)[i1] != other[i2]) {
    if ((*this)[i1] < other[i2]) {
      if (++i1 == this->size())
        return false;
    } else {
      if (++i2 == other.size())
        return false;
    }
  }

  return true;
}

void PermGroup::MovedSet::extend(MovedSet const &other)
{
  MovedSet tmp;

  std::set_union(
    begin(), end(), other.begin(), other.end(), std::back_inserter(tmp));

  *this = tmp;
}

std::vector<PermGroup::EquivalenceClass>
PermGroup::disjoint_decomp_find_equivalence_classes() const
{
  TIMER_START("disjoint decomp find equiv classes");

  std::vector<EquivalenceClass> equivalence_classes;

  MovedSet moved;

  for (Perm const &perm : generators()) {
    moved.init(perm);

    if (equivalence_classes.empty()) {
      equivalence_classes.emplace_back(perm, moved);
      DBG(TRACE) << "Initial equivalence class: {" << perm << "}";
      DBG(TRACE) << "'moved' set is: " << moved;

      continue;
    }

    bool new_class = true;

    for (auto &ec : equivalence_classes) {
      if (!moved.equivalent(ec.moved))
        continue;

      ec.generators.insert(perm);
      DBG(TRACE) << "Updated Equivalence class to " << ec.generators;

      ec.moved.extend(moved);
      DBG(TRACE) << "Updated 'moved' set to " << ec.moved;

      new_class = false;

      break;
    }

    if (new_class) {
      equivalence_classes.emplace_back(perm, moved);
      DBG(TRACE) << "New equivalence class: {" << perm << "}";
      DBG(TRACE) << "'moved' set is: " << moved;
    }
  }

  TIMER_STOP("disjoint decomp find equiv classes");

  return equivalence_classes;
}

void PermGroup::disjoint_decomp_merge_equivalence_classes(
  std::vector<EquivalenceClass> &equivalence_classes) const
{
  TIMER_START("disjoint decomp merge equiv classes");

  unsigned moved_total = 0u;

  for (auto i = 0u; i < equivalence_classes.size(); ++i) {
    EquivalenceClass &ec1 = equivalence_classes[i];

    if (ec1.merged)
      continue;

    for (auto j = i + 1u; j < equivalence_classes.size(); ++j) {
      EquivalenceClass &ec2 = equivalence_classes[j];

      if (ec1.moved.equivalent(ec2.moved)) {
        DBG(TRACE) << "Merging equivalence class " << ec2.generators
                   << " into " << ec1.generators;

        ec1.generators.insert(ec2.generators.begin(), ec2.generators.end());

        ec1.moved.extend(ec2.moved);

        ec2.merged = true;
      }
    }

    if ((moved_total += ec1.moved.size()) == degree())
      break;
  }

  TIMER_STOP("disjoint decomp merge equiv classes");
}

std::vector<PermGroup> PermGroup::disjoint_decomp_incomplete() const
{
  DBG(DEBUG) << "Finding (incomplete) disjoint subgroup decomposition for:";
  DBG(DEBUG) << *this;

  auto equivalence_classes(disjoint_decomp_find_equivalence_classes());

  disjoint_decomp_merge_equivalence_classes(equivalence_classes);

  TIMER_START("disjoint decomp construct groups");

  std::vector<PermGroup> decomp;
  for (auto j = 0u; j < equivalence_classes.size(); ++j) {
    if (equivalence_classes[j].merged)
      continue;

    decomp.emplace_back(degree(), equivalence_classes[j].generators);
  }

  TIMER_STOP("disjoint decomp construct groups");

  DBG(DEBUG) << "Disjoint subgroup generators are:";
#ifndef NDEBUG
  for (PermGroup const &pg : decomp)
    DBG(DEBUG) << pg.generators();
#endif

  return decomp;
}

} // namespace internal

} // namespace mpsym
