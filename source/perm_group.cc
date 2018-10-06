#include <algorithm>
#include <cassert>
#include <climits>
#include <ctime>
#include <random>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

#include "block_system.h"
#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "schreier_sims.h"
#include "util.h"

/**
 * @file perm_group.cc
 * @brief Implements `PermGroup`.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

namespace
{

#ifndef NDEBUG
void debug_orbit_decomposition(std::vector<unsigned> const &orbit_ids,
                                      unsigned n_orbits, unsigned n)
{
  Dbg(Dbg::TRACE) << n_orbits << " orbit(s):";

  std::vector<std::vector<unsigned>> orbits(n_orbits);
  for (unsigned x = 1u; x <= n; ++x) {
    if (orbit_ids[x - 1u])
      orbits[orbit_ids[x - 1u] - 1u].push_back(x);
  }

  for (auto const &orbit : orbits)
    Dbg(Dbg::TRACE) << orbit;
}
#endif

std::vector<PermGroup> disjoint_decomposition_complete_recursive(
  PermGroup const &pg, std::vector<unsigned> const &orbit_ids, unsigned n_orbits)
{
  // iterate over all possible partitions of the set of all orbits into two sets
  assert(n_orbits < 8 * sizeof(unsigned long long));

  for (unsigned long long part = 1ULL; !(part & (1LLU << (n_orbits - 1u))); ++part) {
    std::vector<unsigned> orbit_partition(pg.degree(), 0u);

    for (unsigned x = 0u; x < pg.degree(); ++x) {
      if (!orbit_ids[x])
        continue;

      if ((1ULL << (orbit_ids[x] - 1u)) & part)
        orbit_partition[x] = 2u;
      else
        orbit_partition[x] = 1u;
    }

    Dbg(Dbg::TRACE) << "Considering orbit partition:";

#ifndef NDEBUG
    std::vector<unsigned> moved1, moved2;
    for (unsigned x = 1u; x <= pg.degree(); ++x) {
       if (orbit_partition[x - 1u] == 1u)
         moved1.push_back(x);
       else if (orbit_partition[x - 1u] == 2u)
         moved2.push_back(x);
    }
    Dbg(Dbg::TRACE) << moved1 << " / " << moved2;
#endif

    // determine restricted subgroups
    std::vector<Perm> restricted_gens1;
    std::vector<Perm> restricted_gens2;

    bool recurse = true;
    for (Perm const &gen : pg.bsgs().strong_generators) {
      Dbg(Dbg::TRACE) << "Generator: " << gen;

      // decompose generator into disjoint cycles
      std::vector<unsigned> restricted_gen1(pg.degree());
      std::vector<unsigned> restricted_gen2(pg.degree());
      for (auto i = 1u; i <= pg.degree(); ++i) {
        restricted_gen1[i - 1u] = i;
        restricted_gen2[i - 1u] = i;
      }

      unsigned first = 1u, current = 1u, last = 1u;

      unsigned current_partition = orbit_partition[0u];

      std::vector<bool> processed(pg.degree(), false);

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
            if (next == pg.degree() + 1u) {
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

      Dbg(Dbg::TRACE) << "Restricted generators: "
                      << restricted_gen1 << " / " << restricted_gen2;

      // if generator restricted to orbit partition is not a group member...
      if (!pg.contains_element(Perm(restricted_gen1)) ||
          !pg.contains_element(Perm(restricted_gen2))) {

        recurse = false;

        Dbg(Dbg::TRACE)
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
        Dbg(Dbg::TRACE)
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
      std::vector<unsigned> orbit_ids1(pg.degree(), 0u);
      std::vector<unsigned> orbit_ids2(pg.degree(), 0u);

      std::vector<unsigned> orbit_id_mapping(n_orbits, 0u);
      unsigned current_orbit1 = 0u;
      unsigned current_orbit2 = 0u;

      for (unsigned x = 0u; x < pg.degree(); ++x) {
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

      Dbg(Dbg::TRACE) << "Decomposition is:";
#ifndef NDEBUG
      debug_orbit_decomposition(orbit_ids1, n_orbits1, pg.degree());
      debug_orbit_decomposition(orbit_ids2, n_orbits2, pg.degree());
#endif

      // recurse for both orbit partition elements and return combined result
      std::vector<PermGroup> res1 = disjoint_decomposition_complete_recursive(
        PermGroup(pg.degree(), restricted_gens1), orbit_ids1, n_orbits1);

      std::vector<PermGroup> res2 = disjoint_decomposition_complete_recursive(
        PermGroup(pg.degree(), restricted_gens2), orbit_ids2, n_orbits2);

      for (PermGroup const &perm_group : res2)
        res1.push_back(perm_group);

      return res1;
    }
  }

  Dbg(Dbg::TRACE) << "No further decomposition possible, returning group";
  return {pg};
}

bool orbits_dependent(PermGroup const &pg,
  std::vector<unsigned> orbit1, std::vector<unsigned> orbit2)
{
  Dbg(Dbg::TRACE) << "== Testing whether orbits "
                  << orbit1 << " and " << orbit2 << " are dependent";

  std::unordered_set<Perm> restricted_stabilizers, restricted_elements;

  for (Perm const &perm : pg) {
    Dbg(Dbg::TRACE) << "Looking at generators " << perm;

    bool stabilizes = true;
    for (unsigned o : orbit2) {
      if (perm[o] != o) {
        stabilizes = false;
        break;
      }
    }

    if (stabilizes) {
      Dbg(Dbg::TRACE) << perm << " stabilizes " << orbit2;

      Perm restricted_stabilizer(perm.restricted(orbit1));
      Dbg(Dbg::TRACE) << "Restricted stabilizer is: " << restricted_stabilizer;

      if (!restricted_stabilizer.id())
        restricted_stabilizers.insert(restricted_stabilizer);
    }

    Perm restricted_element(perm.restricted(orbit1));
    Dbg(Dbg::TRACE) << "Restricted group element is: " << restricted_element;

    if (!restricted_element.id())
      restricted_elements.insert(restricted_element);
  }

  bool res = restricted_stabilizers.size() < restricted_elements.size();
  Dbg(Dbg::TRACE) << "=> Orbits " << (res ? "are" : "are not") << " dependent";

  return res;
}

unsigned generate_dependency_classes(
  PermGroup const &pg, std::vector<unsigned> &orbit_ids, unsigned n_orbits)
{
  std::vector<std::vector<unsigned>> orbits(n_orbits);
  for (auto i = 1u; i <= pg.degree(); ++i)
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

      if (orbits_dependent(pg, orbits[i], orbits[j])) {
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

} // anonymous namespace

PermGroup::PermGroup(
  unsigned degree, std::vector<Perm> const &generators,
  ConstructionMethod construction_method,
  TransversalStorageMethod storage_method) : _n(degree)
{
#ifndef NDEBUG
  for (auto const &gen : generators)
    assert(gen.degree() == _n && "all elements have same degree as group");
#endif

  _bsgs.strong_generators = generators;

  if (_bsgs.strong_generators.size() > 0u) {
    switch (construction_method) {
      case SCHREIER_SIMS:
        switch (storage_method) {
          // TODO
          default:
            schreier_sims::schreier_sims<SchreierTree>(_bsgs);
        }
        break;
      case SCHREIER_SIMS_RANDOM:
        switch (storage_method) {
          // TODO
          default:
            schreier_sims::schreier_sims_random<SchreierTree>(_bsgs);
        }
        break;
      default:
        BSGS tmp(BSGS::solve(generators));
        if (!tmp.base.empty()) {
          _bsgs = tmp;
        } else {
          switch (storage_method) {
            // TODO
            default:
              schreier_sims::schreier_sims<SchreierTree>(_bsgs);
          }
        }
    }
  }

  _order = 1u;
  for (unsigned i = 0u; i < _bsgs.base.size(); ++i)
    _order *= _bsgs.orbit(i).size();
}

bool PermGroup::operator==(PermGroup const &rhs) const
{
  assert(_n == rhs.degree() && "comparing permutation groups of equal degree");

  if (_order != rhs.order())
    return false;

  for (Perm const &gen : rhs.bsgs().strong_generators) {
    if (!contains_element(gen))
      return false;
  }

  return true;
}

bool PermGroup::operator!=(PermGroup const &rhs) const
{
  return !(*this == rhs);
}

PermGroup PermGroup::symmetric(unsigned degree)
{
  assert(degree > 0u);

  if (degree == 1u)
    return PermGroup(1u, {Perm(1u)});

  std::vector<unsigned> gen;
  for (unsigned i = 1u; i <= degree; ++i)
    gen.push_back(i);

  return PermGroup(degree, {Perm(degree, {{1, 2}}), Perm(degree, {gen})});
}

PermGroup PermGroup::symmetric(std::vector<unsigned> const &support)
{
  assert(support.size() > 1u);

  unsigned degree = *std::max_element(support.begin(), support.end());

  if (support.size() == 1u)
    return PermGroup(degree, {});

  return PermGroup(degree, {Perm(degree, {{support[0], support[1]}}),
                            Perm(degree, {support})});
}

PermGroup PermGroup::cyclic(unsigned degree)
{
  assert(degree > 0u);

  std::vector<unsigned> gen;
  for (unsigned i = 1u; i <= degree; ++i)
    gen.push_back(i);

  return PermGroup(degree, {Perm(degree, {gen})});
}

PermGroup PermGroup::cyclic(std::vector<unsigned> const &support)
{
  assert(support.size() > 1u);

  unsigned degree = *std::max_element(support.begin(), support.end());

  return PermGroup(degree, {Perm(degree, {support})});
}

PermGroup PermGroup::alternating(unsigned degree)
{
  assert(degree > 0u);

  if (degree == 1u)
    return PermGroup(1, {});

  if (degree == 2u)
    return PermGroup(2, {});

  std::vector<Perm> gens;
  for (unsigned i = 3u; i <= degree; ++i)
    gens.push_back(Perm(degree, {{1, 2, i}}));

  return PermGroup(degree, gens);
}

PermGroup PermGroup::alternating(std::vector<unsigned> const &support)
{
  unsigned degree = *std::max_element(support.begin(), support.end());

  std::vector<Perm> gens;
  for (unsigned i = 2u; i < support.size(); ++i)
    gens.push_back(Perm(degree, {{support[0], support[1], support[i]}}));

  return PermGroup(degree, gens);
}

PermGroup PermGroup::dihedral(unsigned degree)
{
  assert(degree > 0u);

  if (degree == 1u)
    return PermGroup(2, {Perm({2, 1})});

  if (degree == 2u)
    return PermGroup(4, {Perm({2, 1, 3, 4}), Perm({1, 2, 4, 3})});

  std::vector<unsigned> rotation(degree);
  for (unsigned i = 0u; i < degree - 1u; ++i)
    rotation[i] = i + 2u;
  rotation.back() = 1u;

  std::vector<unsigned> reflection(degree);
  for (unsigned i = 0u; i < degree / 2; ++i) {
    reflection[i] = degree - i;
    reflection[degree - i - 1u] = i + 1u;
  }

  if (degree % 2u == 1u)
    reflection[degree / 2u] = degree / 2u + 1u;

  return PermGroup(degree, {Perm(rotation), Perm(reflection)});
}

PermGroup PermGroup::dihedral(std::vector<unsigned> const &support)
{
  assert(support.size() > 2u);

  unsigned degree = *std::max_element(support.begin(), support.end());

  std::vector<unsigned> rotation(degree);
  for (unsigned i = 1u; i <= degree; ++i)
    rotation[i - 1u] = i;

  std::vector<unsigned> reflection(rotation);

  for (unsigned i = 1u; i < support.size(); ++i)
    rotation[support[i - 1u] - 1u] = support[i];
  rotation[support.back() - 1u] = support[0];

  for (unsigned i = 0u; i < support.size() / 2; ++i) {
    reflection[support[i] - 1u] = support[support.size() - i - 1u];
    reflection[support[support.size() - i - 1u] - 1u] = support[i];;
  }

  return PermGroup(degree, {Perm(rotation), Perm(reflection)});
}

PermGroup PermGroup::direct_product(PermGroup const &lhs, PermGroup const &rhs)
{
  auto generators_lhs(std::vector<Perm>(lhs.bsgs().strong_generators));
  auto generators_rhs(std::vector<Perm>(rhs.bsgs().strong_generators));

  unsigned n = lhs.degree() + rhs.degree();

  std::vector<Perm> generators;

  for (auto const &perm : generators_lhs)
    generators.push_back(perm.extended(n));

  for (auto const &perm : generators_rhs)
    generators.push_back(perm.shifted(lhs.degree()));

  return PermGroup(n, generators);
}

bool PermGroup::is_symmetric() const
{
  if (_n == 1u)
    return true;

  unsigned buf = 1u;
  for (unsigned i = _n; i > 0; --i) {
    assert(buf <= UINT_MAX / i);
    buf *= i;
  }

  return _order == buf;
}

bool PermGroup::is_alternating() const
{
  if (_n == 1u)
    return false;

  if (_n == 2u)
    return _order == 1u;

  unsigned buf = 1u;
  for (unsigned i = _n; i > 2; --i) {
    assert(buf <= UINT_MAX / i);
    buf *= i;
  }

  return _order == buf;
}

bool PermGroup::is_transitive() const
{
  std::vector<int> orbit_indices(_n + 1u, -1);
  orbit_indices[1] = 0;

  unsigned processed = 1u;

  for (auto i = 1u; i <= _n; ++i) {
    int orbit_index1 = orbit_indices[i];
    if (orbit_index1 == -1)
      return false;

    for (Perm const &gen : _bsgs.strong_generators) {
      unsigned j = gen[i];

      int orbit_index2 = orbit_indices[j];
      if (orbit_index2 == -1) {
        orbit_indices[j] = orbit_index1;

        if (++processed == _n)
          return true;
      }
    }
  }

  throw std::logic_error("unreachable");
}

std::vector<std::vector<unsigned>> PermGroup::orbits() const
{
  return schreier_sims::orbits(_bsgs.strong_generators);
}

bool PermGroup::contains_element(Perm const &perm) const
{
  assert(perm.degree() == _n && "element has same degree as group");

  auto strip_result = schreier_sims::strip(perm, _bsgs);

  bool ret = (std::get<1>(strip_result) == _bsgs.base.size() + 1) &&
             (std::get<0>(strip_result).id());

  return ret;
}

Perm PermGroup::random_element() const
{
  static std::default_random_engine gen(time(0));

  Perm result(_n);
  for (unsigned i = 0u; i < _bsgs.base.size(); ++i) {
    std::vector<unsigned> orbit = _bsgs.orbit(i);
    std::uniform_int_distribution<> d(0u, orbit.size() - 1u);
    result *= _bsgs.schreier_structures[i]->transversal(orbit[d(gen)]);
  }

  return result;
}

std::vector<PermGroup> PermGroup::disjoint_decomposition(
  bool complete, bool disjoint_orbit_optimization) const
{
  if (complete) {
    return disjoint_decomposition_complete(disjoint_orbit_optimization);
  } else {
#ifndef NDEBUG
    if (disjoint_orbit_optimization) {
      Dbg(Dbg::WARN)
        << "Disjoint orbit optimization ignored during incomplete decomposition";
    }
#endif
    return disjoint_decomposition_incomplete();
  }
}

std::vector<PermGroup> PermGroup::disjoint_decomposition_incomplete() const
{
  Dbg(Dbg::DBG) << "Finding (incomplete) disjoint subgroup decomposition for:";
  Dbg(Dbg::DBG) << *this;

  struct EquivalenceClass {
    std::vector<Perm> generators;
    std::vector<unsigned> moved;
  };

  std::vector<EquivalenceClass> equivalence_classes;

  auto equivalent = [](std::vector<unsigned> const &m1,
                       std::vector<unsigned> const &m2) {

    std::vector<unsigned>::size_type i1 = 0u, i2 = 0u;

    while (m1[i1] != m2[i2]) {
      if (m1[i1] < m2[i2]) {
        if (++i1 == m1.size())
          return false;
      } else {
        if (++i2 == m2.size())
          return false;
      }
    }

    return true;
  };

  auto moved_union = [](std::vector<unsigned> const &m1,
                        std::vector<unsigned> const &m2) {

    std::vector<unsigned> new_moved;

    std::set_union(m1.begin(), m1.end(), m2.begin(), m2.end(),
                   std::back_inserter(new_moved));

    return new_moved;
  };

  std::vector<unsigned> moved;
  for (Perm const &perm : _bsgs.strong_generators) {
    moved.clear();
    for (unsigned i = 1u; i <= perm.degree(); ++i) {
      if (perm[i] != i)
        moved.push_back(i);
    }

    if (equivalence_classes.empty()) {
      equivalence_classes.push_back({{perm}, moved});
      Dbg(Dbg::TRACE) << "Initial equivalence class: " << std::vector<Perm>({perm});
      Dbg(Dbg::TRACE) << "'moved' set is: " << moved;
      continue;
    }

    bool new_class = true;
    for (auto &ec : equivalence_classes) {
      if (!equivalent(moved, ec.moved))
        continue;

      ec.generators.push_back(perm);
      Dbg(Dbg::TRACE) << "Updated Equivalence class to " << ec.generators;

      ec.moved = moved_union(ec.moved, moved);
      Dbg(Dbg::TRACE) << "Updated 'moved' set to " << ec.moved;

      new_class = false;
      break;
    }

    if (new_class) {
      Dbg(Dbg::TRACE) << "New equivalence class containing element: " << perm;
      Dbg(Dbg::TRACE) << "'moved' set is: " << moved;
      equivalence_classes.push_back({{perm}, moved});
    }
  }

  std::vector<bool> merged(equivalence_classes.size(), false);
  unsigned moved_total = 0u;

  std::vector<EquivalenceClass>::size_type i;
  for (i = 0u; i < equivalence_classes.size(); ++i) {
    if (merged[i])
      continue;

    EquivalenceClass &ec1 = equivalence_classes[i];

    for (auto j = i + 1u; j < equivalence_classes.size(); ++j) {
      EquivalenceClass &ec2 = equivalence_classes[j];

      if (equivalent(ec1.moved, ec2.moved)) {
        Dbg(Dbg::TRACE) << "Merging equivalence class " << ec2.generators
                        << " into " << ec1.generators;

        ec1.generators.insert(ec1.generators.end(),
                              ec2.generators.begin(), ec2.generators.end());

        ec1.moved = moved_union(ec1.moved, ec2.moved);

        merged[j] = true;
      }
    }

    if ((moved_total += ec1.moved.size()) == _n)
      break;
  }

  std::vector<PermGroup> decomp;
  for (auto j = 0u; j < equivalence_classes.size(); ++j) {
    if (merged[j])
      continue;

    decomp.push_back(PermGroup(_n, equivalence_classes[j].generators));
  }

  Dbg(Dbg::DBG) << "Disjunct subgroup generators are:";
#ifndef NDEBUG
  for (PermGroup const &pg : decomp)
    Dbg(Dbg::DBG) << pg.bsgs().strong_generators;
#endif

  return decomp;
}

std::vector<PermGroup> PermGroup::disjoint_decomposition_complete(
  bool disjoint_orbit_optimization) const
{
  Dbg(Dbg::DBG) << "Finding (complete) disjoint subgroup decomposition for:";
  Dbg(Dbg::DBG) << *this;

  // determine complete orbit decomposition
  std::vector<unsigned> orbit_ids(_n, 0u);
  orbit_ids[0u] = 1u;

  unsigned n_processed = 1u;
  unsigned n_orbits = 1u;
  for (unsigned x = 1u; x <= _n; ++x) {
    unsigned orbit_id = 0u;
    std::vector<unsigned> new_orbit_elems;

    if (!orbit_ids[x - 1u])
      new_orbit_elems.push_back(x);
    else
      orbit_id = orbit_ids[x - 1u];

    for (Perm const &gen : _bsgs.strong_generators) {
      unsigned y = gen[x];

      if (y != x) {
        if (!orbit_ids[y - 1u]) {
          auto it =
            std::find(new_orbit_elems.begin(), new_orbit_elems.end(), y);

          if (it == new_orbit_elems.end())
            new_orbit_elems.push_back(y);
        } else if (!orbit_id)
          orbit_id = orbit_ids[y - 1u];
      }
    }

    if (!orbit_id)
      orbit_id = ++n_orbits;

    bool done = false;
    for (unsigned y : new_orbit_elems) {
      orbit_ids[y - 1u] = orbit_id;

      if (++n_processed == _n) {
        done = true;
        break;
      }
    }

    if (done)
      break;
  }

  Dbg(Dbg::TRACE) << "=== Orbit decomposition:";
#ifndef NDEBUG
  debug_orbit_decomposition(orbit_ids, n_orbits, _n);
#endif

  if (disjoint_orbit_optimization) {
    Dbg(Dbg::TRACE) << "=== Using dependent orbit optimization";
    n_orbits = generate_dependency_classes(*this, orbit_ids, n_orbits);

    Dbg(Dbg::TRACE) << "==> Grouped dependency class unions:";
#ifndef NDEBUG
    debug_orbit_decomposition(orbit_ids, n_orbits, _n);
#endif
  }

  std::vector<PermGroup> decomp =
    disjoint_decomposition_complete_recursive(*this, orbit_ids, n_orbits);

  Dbg(Dbg::DBG) << "=== Disjoint subgroup generators";
  for (PermGroup const &pg : decomp)
    Dbg(Dbg::DBG) << pg.bsgs().strong_generators;

  return decomp;
}

std::vector<PermGroup> PermGroup::wreath_decomposition() const
{
  Dbg(Dbg::DBG) << "Finding wreath product decomposition for:";
  Dbg(Dbg::DBG) << *this;

  auto blocksystems(BlockSystem::non_trivial(*this));
  auto gens(_bsgs.strong_generators);

  for (BlockSystem const &bs : blocksystems) {
    Dbg(Dbg::TRACE) << "Considering block system: " << bs;
    unsigned d = bs.size();

    PermGroup block_permuter = bs.block_permuter(gens);
    Dbg(Dbg::TRACE) << "Block permuter is:";
    Dbg(Dbg::TRACE) << block_permuter;

    std::vector<std::vector<Perm>> sigma(d);

    sigma[0] = BlockSystem::block_stabilizers(gens, bs[0]);
    Dbg(Dbg::TRACE) << "Block stabilizer of " << bs[0] << " is: " << sigma[0];

    unsigned tmp = PermGroup(_n, sigma[0]).order();
    if (_order != pow(tmp, d) * block_permuter.order()) {
      Dbg(Dbg::TRACE)
        << "Group order equality not satisfied, skipping block system";
      continue;
    }

    Dbg(Dbg::TRACE) << "Stabilizers of remaining blocks are:";
    for (unsigned i = 1u; i < d; ++i) {
      sigma[i] = BlockSystem::block_stabilizers(gens, bs[i]);
      Dbg(Dbg::TRACE) << bs[i] << " => " << sigma[i];
    }

    // try to find monomorphism from block permuter to group using heuristic
    std::vector<Perm> block_permuter_generator_image;

    std::vector<unsigned> tmp_perm(_n);
    for (Perm const &gen : block_permuter.bsgs().strong_generators) {
      for (unsigned i = 0u; i < d; ++i) {
        auto block(bs[i]);

        for (auto j = 0u; j < block.size(); ++j)
          tmp_perm[block[j] - 1u] = bs[gen[i + 1u] - 1u][j];
      }

      block_permuter_generator_image.push_back(Perm(tmp_perm));
    }

    Dbg(Dbg::TRACE) << "Heuristic monomorphism image is:";
    Dbg(Dbg::TRACE) << block_permuter_generator_image;

    bool found_monomorphism = true;

    auto classes(bs.classes());

    std::vector<Perm> block_permuter_generator_reconstruction;

    tmp_perm.resize(d);
    for (Perm const &gen : block_permuter_generator_image) {
      for (unsigned i = 0u; i < d; ++i) {
        unsigned x = bs[i][0];
        unsigned y = gen[x];

        tmp_perm[i] = classes[y - 1u];
      }

      Perm reconstructed_gen(tmp_perm);

      if (!block_permuter.contains_element(reconstructed_gen)) {
        found_monomorphism = false;
        break;
      }

      block_permuter_generator_reconstruction.push_back(reconstructed_gen);
    }

    Dbg(Dbg::TRACE) << "Block permuter reconstruction yields generators:";
    Dbg(Dbg::TRACE) << block_permuter_generator_reconstruction;

    if (found_monomorphism) {
      PermGroup block_permuter_reconstructed(
        d, block_permuter_generator_reconstruction);

      if (block_permuter_reconstructed.order() != block_permuter.order())
        found_monomorphism = false;
    }

    if (!found_monomorphism) {
      Dbg(Dbg::WARN)
        << "Wreath decomposition exists but was not found by heuristic";

      break;
    }

    std::vector<PermGroup> res(d + 1u);

    res[0] = PermGroup(_n, block_permuter_generator_image);
    for (unsigned i = 0u; i < d; ++i)
      res[i + 1u] = PermGroup(_n, sigma[i]);

    Dbg(Dbg::TRACE)
      << "==> Found wreath product decomposition, listing generators:";
#ifndef NDEBUG
    for (PermGroup const &pg : res)
      Dbg(Dbg::TRACE) << pg.bsgs().strong_generators;
#endif

    return res;
  }

  Dbg(Dbg::TRACE) << "==> No wreath product decomposition found";
  return std::vector<PermGroup>();
}

PermGroup::const_iterator::const_iterator(PermGroup const &pg)
  : _trivial(pg.bsgs().base.empty()), _end(false)
{
  if (_trivial) {
    _current_result = Perm(pg.degree());
  } else {
    for (unsigned i = 0u; i < pg.bsgs().base.size(); ++i) {
      _state.push_back(0u);
      _transversals.push_back(pg.bsgs().transversals(i));
      _current_factors.push_back(_transversals.back()[0]);
    }

    update_result();
  }
}

PermGroup::const_iterator PermGroup::const_iterator::operator++()
{
  PermGroup::const_iterator pre(*this);
  next_state();
  return pre;
}

bool PermGroup::const_iterator::operator==(
  PermGroup::const_iterator const &rhs) const
{
  if (_end != rhs._end)
    return false;

  if (_end && rhs._end)
    return true;

  for (unsigned i = 0u; i < _state.size(); ++i) {
    if (_state[i] != rhs._state[i])
      return false;
  }
  return true;
}

void PermGroup::const_iterator::next_state()
{
  if (_trivial) {
    _end = true;
    return;
  }

  for (unsigned i = 0u; i < _state.size(); ++i) {
    _state[i]++;
    if (_state[i] == _transversals[i].size())
      _state[i] = 0u;

    _current_factors[i] = _transversals[i][_state[i]];

    if (i == _state.size() - 1u && _state[i] == 0u) {
      _end = true;
      break;
    }

    if (_state[i] != 0u)
      break;
  }

  update_result();
}

void PermGroup::const_iterator::update_result()
{
  _current_result = _current_factors[0];
  for (unsigned j = 1u; j < _current_factors.size(); ++j)
    _current_result = _current_factors[j] * _current_result;
}

std::ostream& operator<<(std::ostream& stream, PermGroup const &pg) {
  stream << pg.bsgs();
  return stream;
}

} // namespace cgtl
