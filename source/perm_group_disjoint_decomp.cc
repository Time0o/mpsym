#include <algorithm>
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

#ifndef NDEBUG
namespace
{

void debug_orbit_decomposition(std::vector<unsigned> const &orbit_ids,
                               unsigned n_orbits,
                               unsigned n)
{
  DBG(TRACE) << n_orbits << " orbit(s):";

  std::vector<std::vector<unsigned>> orbits(n_orbits);
  for (unsigned x = 1u; x <= n; ++x) {
    if (orbit_ids[x - 1u])
      orbits[orbit_ids[x - 1u] - 1u].push_back(x);
  }

  for (auto const &orbit : orbits)
    DBG(TRACE) << orbit;
}

} // namespace
#endif

namespace cgtl
{

std::vector<PermGroup> PermGroup::disjoint_decomposition(
  bool complete,
  bool disjoint_orbit_optimization) const
{
  if (complete) {
    return disjoint_decomp_complete(disjoint_orbit_optimization);
  } else {
#ifndef NDEBUG
    if (disjoint_orbit_optimization) {
      DBG(WARN)
        << "Disjoint orbit optimization ignored during incomplete decomposition";
    }
#endif
    return disjoint_decomp_incomplete();
  }
}

bool PermGroup::disjoint_decomp_complete_orbits_dependent(
  std::vector<unsigned> orbit1,
  std::vector<unsigned> orbit2) const
{
  DBG(TRACE) << "== Testing whether orbits "
             << orbit1 << " and " << orbit2 << " are dependent";

  std::unordered_set<Perm> restricted_stabilizers, restricted_elements;

  for (Perm const &perm : *this) {
    DBG(TRACE) << "Looking at generators " << perm;

    bool stabilizes = true;
    for (unsigned o : orbit2) {
      if (perm[o] != o) {
        stabilizes = false;
        break;
      }
    }

    if (stabilizes) {
      DBG(TRACE) << perm << " stabilizes " << orbit2;

      Perm restricted_stabilizer(perm.restricted(orbit1));
      DBG(TRACE) << "Restricted stabilizer is: " << restricted_stabilizer;

      if (!restricted_stabilizer.id())
        restricted_stabilizers.insert(restricted_stabilizer);
    }

    Perm restricted_element(perm.restricted(orbit1));
    DBG(TRACE) << "Restricted group element is: " << restricted_element;

    if (!restricted_element.id())
      restricted_elements.insert(restricted_element);
  }

  bool res = restricted_stabilizers.size() < restricted_elements.size();
  DBG(TRACE) << "=> Orbits " << (res ? "are" : "are not") << " dependent";

  return res;
}


unsigned PermGroup::disjoint_decomp_complete_generate_dependency_classes(
  std::vector<unsigned> &orbit_ids,
  unsigned n_orbits) const
{
  std::vector<std::vector<unsigned>> orbits(n_orbits);
  for (auto i = 1u; i <= degree(); ++i)
    orbits[orbit_ids[i - 1u] - 1u].push_back(i);

  std::vector<std::vector<unsigned>> dependency_classes;
  unsigned n_dependency_classes = 0u;

  std::vector<bool> processed(n_orbits, false);
  unsigned n_processed = 0u;

  for (auto i = 0u; i < n_orbits; ++i) {
    if (processed[i])
      continue;

    std::vector<unsigned> merge {i + 1u};

    for (auto j = i + 1u; j < n_orbits; ++j) {
      if (processed[j])
        continue;

      if (disjoint_decomp_complete_orbits_dependent(orbits[i], orbits[j])) {
        merge.push_back(j + 1u);
        processed[j] = true;
        ++n_processed;
      }
    }

    ++n_dependency_classes;
    for (unsigned &orbit_id : orbit_ids) {
      if (std::find(merge.begin(), merge.end(), (orbit_id)) != merge.end())
        orbit_id = n_dependency_classes;
    }

    processed[i] = true;
    if (++n_processed == n_orbits)
      return n_dependency_classes;
  }

  return n_dependency_classes;
}

std::vector<PermGroup> PermGroup::disjoint_decomp_complete_recursive(
  std::vector<unsigned> const &orbit_ids,
  unsigned n_orbits,
  PermGroup const *pg) const
{
  if (!pg)
    pg = this;

  // iterate over all possible partitions of the set of all orbits into two sets
  assert(n_orbits < 8 * sizeof(unsigned long long));

  for (unsigned long long part = 1ULL; !(part & (1LLU << (n_orbits - 1u))); ++part) {
    std::vector<unsigned> orbit_partition(pg->degree(), 0u);

    for (unsigned x = 0u; x < pg->degree(); ++x) {
      if (!orbit_ids[x])
        continue;

      if ((1ULL << (orbit_ids[x] - 1u)) & part)
        orbit_partition[x] = 2u;
      else
        orbit_partition[x] = 1u;
    }

    DBG(TRACE) << "Considering orbit partition:";

#ifndef NDEBUG
    std::vector<unsigned> moved1, moved2;
    for (unsigned x = 1u; x <= pg->degree(); ++x) {
       if (orbit_partition[x - 1u] == 1u)
         moved1.push_back(x);
       else if (orbit_partition[x - 1u] == 2u)
         moved2.push_back(x);
    }
    DBG(TRACE) << moved1 << " / " << moved2;
#endif

    // determine restricted subgroups
    std::vector<Perm> restricted_gens1;
    std::vector<Perm> restricted_gens2;

    bool recurse = true;
    for (Perm const &gen : pg->generators()) {
      DBG(TRACE) << "Generator: " << gen;

      // decompose generator into disjoint cycles
      std::vector<unsigned> restricted_gen1(pg->degree());
      std::vector<unsigned> restricted_gen2(pg->degree());
      for (auto i = 1u; i <= pg->degree(); ++i) {
        restricted_gen1[i - 1u] = i;
        restricted_gen2[i - 1u] = i;
      }

      unsigned first = 1u, current = 1u, last = 1u;

      unsigned current_partition = orbit_partition[0u];

      std::vector<bool> processed(pg->degree(), false);

      for (;;) {
        processed[current - 1u] = true;
        last = current;
        current = gen[current];

        if (current_partition == 1u)
          restricted_gen1[last - 1u] = current;
        else if (current_partition == 2u)
          restricted_gen2[last - 1u] = current;

        if (current == first) {
          bool done = false;

          unsigned next = first + 1u;
          for (;;) {
            if (next == pg->degree() + 1u) {
              done = true;
              break;
            }
            if (!processed[next - 1u])
              break;

            ++next;
          }

          if (done)
            break;
          else {
            first = next;
            current = first;
            current_partition = orbit_partition[current - 1u];
          }
        }
      }

      DBG(TRACE) << "Restricted generators: "
                 << restricted_gen1 << " / " << restricted_gen2;

      // if generator restricted to orbit partition is not a group member...
      if (!pg->contains_element(Perm(restricted_gen1)) ||
          !pg->contains_element(Perm(restricted_gen2))) {

        recurse = false;

        DBG(TRACE)
          << "Restricted groups are not a disjoint subgroup decomposition";

        break;
      }
      else {
        restricted_gens1.push_back(Perm(restricted_gen1));
        restricted_gens2.push_back(Perm(restricted_gen2));
      }
    }

    // recurse for group restricted to orbit partitions
    if (recurse) {
        DBG(TRACE)
          << "Restricted groups are a disjoint subgroup decomposition, recursing...";

      // compute sizes orbit partition elements
      unsigned n_orbits1 = 0u;
      unsigned n_orbits2 = 0u;

      for (unsigned shift = 0u; shift < n_orbits; ++shift) {
        if (!(part & (1ULL << shift)))
          ++n_orbits1;
        else
          ++n_orbits2;
      }

      // compute orbit partition elements
      std::vector<unsigned> orbit_ids1(pg->degree(), 0u);
      std::vector<unsigned> orbit_ids2(pg->degree(), 0u);

      std::vector<unsigned> orbit_id_mapping(n_orbits, 0u);
      unsigned current_orbit1 = 0u;
      unsigned current_orbit2 = 0u;

      for (unsigned x = 0u; x < pg->degree(); ++x) {
        if (orbit_partition[x] == 1u) {
          if (!orbit_id_mapping[orbit_ids[x] - 1u])
            orbit_id_mapping[orbit_ids[x] - 1u] = ++current_orbit1;

          orbit_ids1[x] = orbit_id_mapping[orbit_ids[x] - 1u];
        }
        else if (orbit_partition[x] == 2u) {
          if (!orbit_id_mapping[orbit_ids[x] - 1u])
            orbit_id_mapping[orbit_ids[x] - 1u] = ++current_orbit2;

          orbit_ids2[x] = orbit_id_mapping[orbit_ids[x] - 1u];
        }
      }

      DBG(TRACE) << "Decomposition is:";
#ifndef NDEBUG
      debug_orbit_decomposition(orbit_ids1, n_orbits1, pg->degree());
      debug_orbit_decomposition(orbit_ids2, n_orbits2, pg->degree());
#endif

      // recurse for both orbit partition elements and return combined result
      PermGroup pg1(pg->degree(),
                    PermSet(restricted_gens1.begin(), restricted_gens1.end()));

      std::vector<PermGroup> res1 = disjoint_decomp_complete_recursive(
        orbit_ids1, n_orbits1, &pg1);

      PermGroup pg2(pg->degree(),
                    PermSet(restricted_gens2.begin(), restricted_gens2.end()));

      std::vector<PermGroup> res2 = disjoint_decomp_complete_recursive(
        orbit_ids2, n_orbits2, &pg2);

      for (PermGroup const &perm_group : res2)
        res1.push_back(perm_group);

      return res1;
    }
  }

  DBG(TRACE) << "No further decomposition possible, returning group";
  return {*pg};
}

std::vector<PermGroup> PermGroup::disjoint_decomp_complete(
  bool disjoint_orbit_optimization) const
{
  DBG(DEBUG) << "Finding (complete) disjoint subgroup decomposition for:";
  DBG(DEBUG) << *this;

  auto tmp(orbit_partition(generators()));
  auto orbit_ids(tmp.first);
  auto n_orbits = tmp.second;

  DBG(TRACE) << "=== Orbit decomposition:";
#ifndef NDEBUG
  debug_orbit_decomposition(orbit_ids, n_orbits, degree());
#endif

  if (disjoint_orbit_optimization) {
    DBG(TRACE) << "=== Using dependent orbit optimization";
    n_orbits = disjoint_decomp_complete_generate_dependency_classes(orbit_ids,
                                                                    n_orbits);

    DBG(TRACE) << "==> Grouped dependency class unions:";
#ifndef NDEBUG
    debug_orbit_decomposition(orbit_ids, n_orbits, degree());
#endif
  }

  std::vector<PermGroup> decomp = disjoint_decomp_complete_recursive(orbit_ids,
                                                                     n_orbits);

  DBG(DEBUG) << "=== Disjoint subgroup generators";
  for (PermGroup const &pg : decomp)
    DBG(DEBUG) << pg.generators();

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
PermGroup::disjoint_decomp_incomplete_find_equivalence_classes() const
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

void PermGroup::disjoint_decomp_incomplete_merge_equivalence_classes(
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

  auto equivalence_classes(
    disjoint_decomp_incomplete_find_equivalence_classes());

  disjoint_decomp_incomplete_merge_equivalence_classes(equivalence_classes);

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

} // namespace cgtl
