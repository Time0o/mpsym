#include <algorithm>
#include <cassert>
#include <climits>
#include <ctime>
#include <random>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "orbits.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"
#include "schreier_structure.h"
#include "util.h"

/**
 * @file perm_group.cc
 * @brief Implements `PermGroup`.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

PermGroup::PermGroup(unsigned degree,
                     PermSet const &generators,
                     BSGS::Construction construction,
                     BSGS::Transversals transversals)
{
  if (generators.empty() || (generators.size() == 1u && generators[0].id())) {
    _bsgs = BSGS(degree);

    _order = 1ULL;

  } else {
    _bsgs = BSGS(degree, generators, construction, transversals);

    _order = 1ULL;
    for (unsigned i = 0u; i < _bsgs.base_size(); ++i) {
      unsigned long long orbit_size = _bsgs.orbit(i).size();

      // TODO: this restriction might need to be lifted
      assert(_order <= ULLONG_MAX / orbit_size &&
             "group order representable by unsigned long long");

      _order *= orbit_size;
    }
  }
}

bool PermGroup::operator==(PermGroup const &rhs) const
{
  assert(rhs.degree() == degree()
         && "comparing permutation groups of equal degree");

  if (_order != rhs.order())
    return false;

  for (Perm const &gen : rhs.generators()) {
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
  // TODO: explicit BSGS

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
  // TODO: explicit BSGS

  assert(support.size() > 1u);

  unsigned degree = *std::max_element(support.begin(), support.end());

  if (support.size() == 1u)
    return PermGroup(degree, {});

  return PermGroup(degree, {Perm(degree, {{support[0], support[1]}}),
                            Perm(degree, {support})});
}

PermGroup PermGroup::cyclic(unsigned degree)
{
  // TODO: explicit BSGS

  assert(degree > 0u);

  std::vector<unsigned> gen;
  for (unsigned i = 1u; i <= degree; ++i)
    gen.push_back(i);

  return PermGroup(degree, {Perm(degree, {gen})});
}

PermGroup PermGroup::cyclic(std::vector<unsigned> const &support)
{
  // TODO: explicit BSGS

  assert(support.size() > 1u);

  unsigned degree = *std::max_element(support.begin(), support.end());

  return PermGroup(degree, {Perm(degree, {support})});
}

PermGroup PermGroup::alternating(unsigned degree)
{
  // TODO: explicit BSGS

  assert(degree > 0u);

  if (degree == 1u)
    return PermGroup(1, {});

  if (degree == 2u)
    return PermGroup(2, {});

  PermSet gens;
  for (unsigned i = 3u; i <= degree; ++i)
    gens.emplace(Perm(degree, {{1, 2, i}}));

  return PermGroup(degree, gens);
}

PermGroup PermGroup::alternating(std::vector<unsigned> const &support)
{
  // TODO: explicit BSGS

  unsigned degree = *std::max_element(support.begin(), support.end());

  PermSet gens;
  for (unsigned i = 2u; i < support.size(); ++i)
    gens.emplace(Perm(degree, {{support[0], support[1], support[i]}}));

  return PermGroup(degree, gens);
}

PermGroup PermGroup::dihedral(unsigned degree)
{
  // TODO: explicit BSGS

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
  // TODO: explicit BSGS

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

PermGroup PermGroup::wreath_product(PermGroup const &lhs_, PermGroup const &rhs_)
{
  auto lhs(lhs_.generators());
  auto rhs(rhs_.generators());

  lhs.assert_not_empty();
  rhs.assert_not_empty();

  lhs.minimize_degree();
  rhs.minimize_degree();

  unsigned degree = lhs.degree() * rhs.degree();

  PermSet wreath_product_generators;

  for (unsigned i = 0u; i < rhs.degree(); ++i) {
    for (Perm const &perm : lhs)
      wreath_product_generators.insert(
        perm.shifted(lhs.degree() * i).extended(degree));
  }

  for (Perm const &gen_rhs : rhs) {
    std::vector<std::vector<unsigned>> cycles {gen_rhs.cycles()};
    for (auto &cycle : cycles) {
      for (unsigned &x : cycle)
        x = (x - 1u) * lhs.degree() + 1u;
    }

    std::vector<std::vector<unsigned>> shifted_cycles {cycles};

    for (unsigned i = 1u; i < lhs.degree(); ++i) {
      for (auto const &cycle : cycles) {
        std::vector<unsigned> shifted_cycle(cycle);

        for (unsigned &x : shifted_cycle)
          x += i;

        shifted_cycles.push_back(shifted_cycle);
      }
    }

    wreath_product_generators.emplace(degree, shifted_cycles);
  }

  return PermGroup(degree, wreath_product_generators);
}

bool PermGroup::is_symmetric() const
{
  if (degree() == 1u)
    return true;

  unsigned buf = 1u;
  for (unsigned i = degree(); i > 0; --i) {
    assert(buf <= UINT_MAX / i);
    buf *= i;
  }

  return _order == buf;
}

bool PermGroup::is_alternating() const
{
  if (degree() == 1u)
    return false;

  if (degree() == 2u)
    return _order == 1u;

  unsigned buf = 1u;
  for (unsigned i = degree(); i > 2; --i) {
    assert(buf <= UINT_MAX / i);
    buf *= i;
  }

  return _order == buf;
}

bool PermGroup::is_transitive() const
{ return Orbit::generate(1u, generators()).size() == degree(); }

bool PermGroup::contains_element(Perm const &perm) const
{
  assert(perm.degree() == degree() && "element has same degree as group");

  return _bsgs.strips_completely(perm);
}

Perm PermGroup::random_element() const
{
  static auto re(util::random_engine());

  Perm result(degree());
  for (unsigned i = 0u; i < _bsgs.base_size(); ++i) {
    std::vector<unsigned> orbit = _bsgs.orbit(i);
    std::uniform_int_distribution<> d(0u, orbit.size() - 1u);
    result *= _bsgs.transversal(i, orbit[d(re)]);
  }

  return result;
}

PermGroup::const_iterator::const_iterator(PermGroup const &pg)
  : _trivial(pg.bsgs().base_empty()), _end(false)
{
  if (_trivial) {
    _current_result = Perm(pg.degree());
  } else {
    for (unsigned i = 0u; i < pg.bsgs().base_size(); ++i) {
      _state.push_back(0u);

      auto transv = pg.bsgs().transversals(i);

      _transversals.push_back(transv);
      _current_factors.insert(transv[0]);
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

std::ostream &operator<<(std::ostream &os, PermGroup const &pg)
{
  os << pg.bsgs() << "\n"
     << "ORDER: " << pg._order;

  return os;
}

} // namespace cgtl
