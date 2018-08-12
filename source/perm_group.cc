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
  Dbg(Dbg::DBG) << "Finding disjoint subgroup decomposition for:";
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

  Dbg(Dbg::DBG) << "Disjunct subgroup generators are:";
  std::vector<PermGroup> decomp;
  for (auto j = 0u; j < equivalence_classes.size(); ++j) {
    if (merged[j])
      continue;

    decomp.push_back(PermGroup(_n, equivalence_classes[j].generators));
    Dbg(Dbg::DBG) << equivalence_classes[j].generators;
  }
  return decomp;
}

std::vector<PermGroup> PermGroup::disjoint_decomposition_complete() const
{
  throw std::runtime_error("not implemented");
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
