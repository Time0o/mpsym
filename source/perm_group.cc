#include <algorithm>
#include <cassert>
#include <ctime>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "schreier_sims.h"

namespace cgtl
{

PermGroup::PermGroup(unsigned degree, std::vector<Perm> const &generators,
  SchreierSims::Variant schreier_sims_method)
  : _n(degree), _bsgs(generators, schreier_sims_method)
{
#ifndef NDEBUG
  for (auto const &gen : generators)
    assert(gen.degree() == _n && "all elements have same degree as group");
#endif

  _order = 1u;
  for (auto const &b : _bsgs)
    _order *= b.orbit().size();
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

PermGroup PermGroup::cyclic(unsigned degree)
{
  assert(degree > 0u);

  std::vector<unsigned> gen;
  for (unsigned i = 1u; i <= degree; ++i)
    gen.push_back(i);

  return PermGroup(degree, {Perm(degree, {gen})});
}

PermGroup PermGroup::alternating(unsigned degree)
{
  assert(degree > 2u);

  std::vector<Perm> gens;
  for (unsigned i = 3u; i <= degree; ++i)
    gens.push_back(Perm(degree, {{1, 2, i}}));

  return PermGroup(degree, gens);
}

bool PermGroup::is_element(Perm const &perm) const
{
  assert(perm.degree() == _n && "element has same degree as group");

  Dbg(Dbg::DBG) << "Performing membership test for " << perm << " in:";
  Dbg(Dbg::DBG) << (*this);

  auto strip_result = _bsgs.strip(perm);

  Dbg(Dbg::TRACE) << "Strip returned " << std::get<0>(strip_result) << ", "
                  << std::get<1>(strip_result);

  bool ret = (std::get<1>(strip_result) == _bsgs.size() + 1) &&
             (std::get<0>(strip_result).id());

  Dbg(Dbg::DBG) << (ret ? "=> Member" : "=> No Member");

  return ret;
}

Perm PermGroup::random_element() const
{
  static std::default_random_engine gen(time(0));

  Perm result(_n);
  for (auto const &b : _bsgs) {
    std::vector<unsigned> orbit = b.orbit();
    std::uniform_int_distribution<> d(0u, orbit.size() - 1u);
    result *= b.transversal(orbit[d(gen)]);
  }

  return result;
}

std::vector<PermGroup> PermGroup::disjoint_decomposition(bool complete) const
{
  if (complete)
    return disjoint_decomposition_complete();
  else
    return disjoint_decomposition_incomplete();
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
  for (Perm const &perm : _bsgs.sgs()) {
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
    Dbg(Dbg::DBG) << pg.bsgs().sgs();
#endif

  return decomp;
}

#ifndef NDEBUG
static void debug_orbit_decomposition(std::vector<unsigned> const &orbit_ids,
                                      unsigned n_orbits, unsigned n)
{
  Dbg(Dbg::TRACE) << n_orbits << " orbits:";

  std::vector<std::vector<unsigned>> orbits(n_orbits);
  for (unsigned x = 1u; x <= n; ++x) {
    if (orbit_ids[x - 1u])
      orbits[orbit_ids[x - 1u] - 1u].push_back(x);
  }

  for (auto const &orbit : orbits)
    Dbg(Dbg::TRACE) << orbit;
}
#endif

static std::vector<PermGroup> disjoint_decomposition_complete_recursive(
  PermGroup const &pg, std::vector<unsigned> const &orbit_ids, unsigned n_orbits)
{
  Dbg(Dbg::TRACE) << "Orbit decomposition:";
#ifndef NDEBUG
  debug_orbit_decomposition(orbit_ids, n_orbits, pg.degree());
#endif

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
    for (Perm const &gen : pg.bsgs().sgs()) {
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
      if (!pg.is_element(Perm(restricted_gen1)) ||
          !pg.is_element(Perm(restricted_gen2))) {

        recurse = false;

        Dbg(Dbg::TRACE)
          << "Restricted groups are not a disjunct subgroup decomposition";

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
          << "Restricted groups are a disjunct subgroup decomposition, recursing...";

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

std::vector<PermGroup> PermGroup::disjoint_decomposition_complete() const
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

    for (Perm const &gen : _bsgs.sgs()) {
      unsigned y = gen[x];

      if (y != x) {
        if (!orbit_ids[y - 1u])
          new_orbit_elems.push_back(y);
        else if (!orbit_id)
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

  std::vector<PermGroup> decomp = disjoint_decomposition_complete_recursive(
    *this, orbit_ids, n_orbits);

  Dbg(Dbg::DBG) << "Disjunct subgroup generators are:";
  for (PermGroup const &pg : decomp)
    Dbg(Dbg::DBG) << pg.bsgs().sgs();

  return decomp;
}

PermGroup::const_iterator::const_iterator(PermGroup const &pg)
  : _trivial(pg.bsgs().trivial()), _end(false)
{
  if (_trivial) {
    _current_result = Perm(pg.degree());
  } else {
    for (auto const &b : pg.bsgs()) {
      _state.push_back(0u);
      _transversals.push_back(b.transversals());
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

} // namespace cgtl
