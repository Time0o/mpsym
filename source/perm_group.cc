#include <algorithm>
#include <cassert>
#include <climits>
#include <ctime>
#include <queue>
#include <random>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
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

namespace
{

unsigned compress_generators(unsigned degree, PermSet &generators)
{
  std::vector<std::vector<unsigned>> moved_sets(generators.size());

  std::vector<unsigned> compression_mapping(degree + 1u);
  for (unsigned i = 1u; i <= degree; ++i)
    compression_mapping[i] = i;

  std::queue<unsigned> non_moved_queue;

  unsigned new_degree = 1u;

  for (unsigned i = 1u; i <= degree; ++i) {
    bool moved = false;
    for (auto j = 0u; j < generators.size(); ++j) {
      if (generators[j][i] != i) {
        moved_sets[j].push_back(i);
        moved = true;
      }
    }

    if (moved) {
      if (!non_moved_queue.empty()) {
        unsigned compress_to = non_moved_queue.front();
        compression_mapping[i] = compress_to;
        new_degree = compress_to;

        non_moved_queue.pop();
        non_moved_queue.push(i);
      } else {
        new_degree = i;
      }
    } else {
      non_moved_queue.push(i);
    }
  }

  std::vector<unsigned> id(new_degree);
  for (unsigned i = 1u; i <= new_degree; ++i)
    id[i - 1u] = i;

  for (unsigned i = 0u; i < generators.size(); ++i) {
    auto gen(id);
    for (unsigned j = 0u; j < moved_sets[i].size(); ++j) {
      unsigned x = moved_sets[i][j];
      unsigned y = generators[i][x];
      gen[compression_mapping[x] - 1u] = compression_mapping[y];
    }

    generators[i] = Perm(gen);
  }

  return new_degree;
}

} // anonymous namespace

PermGroup::PermGroup(unsigned degree,
                     PermSet const &generators,
                     BSGS::Construction construction,
                     BSGS::Transversals transversals)
: _bsgs(degree, generators, construction, transversals)
{
  _order = 1ULL;
  for (unsigned i = 0u; i < _bsgs.base_size(); ++i) {
    unsigned long long orbit_size = _bsgs.orbit(i).size();

    // TODO: this restriction might need to be lifted
    assert(_order <= ULLONG_MAX / orbit_size &&
           "group order representable by unsigned long long");

    _order *= orbit_size;
  }
}

bool PermGroup::operator==(PermGroup const &rhs) const
{
  assert(rhs.degree() == degree()
         && "comparing permutation groups of equal degree");

  if (_order != rhs.order())
    return false;

  for (Perm const &gen : rhs.bsgs().strong_generators()) {
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

  PermSet gens;
  for (unsigned i = 3u; i <= degree; ++i)
    gens.emplace(Perm(degree, {{1, 2, i}}));

  return PermGroup(degree, gens);
}

PermGroup PermGroup::alternating(std::vector<unsigned> const &support)
{
  unsigned degree = *std::max_element(support.begin(), support.end());

  PermSet gens;
  for (unsigned i = 2u; i < support.size(); ++i)
    gens.emplace(Perm(degree, {{support[0], support[1], support[i]}}));

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
  return direct_product(lhs.bsgs().strong_generators(),
                        rhs.bsgs().strong_generators());
}

PermGroup PermGroup::direct_product(PermSet const &lhs, PermSet const &rhs)
{
  lhs.assert_not_empty();
  rhs.assert_not_empty();

  unsigned degree = lhs.degree() + rhs.degree();

  PermSet generators;

  for (auto const &perm : lhs)
    generators.insert(perm.extended(degree));

  for (auto const &perm : rhs)
    generators.insert(perm.shifted(lhs.degree()));

  return PermGroup(degree, generators);
}

PermGroup PermGroup::wreath_product(PermGroup const &lhs, PermGroup const &rhs)
{
  return wreath_product(lhs.bsgs().strong_generators(),
                        rhs.bsgs().strong_generators());
}

PermGroup PermGroup::wreath_product(PermSet const &lhs, PermSet const &rhs)
{
  lhs.assert_not_empty();
  rhs.assert_not_empty();

  unsigned degree_lhs = lhs.degree();
  unsigned degree_rhs = rhs.degree();

  auto lhs_compressed(lhs);
  degree_lhs = compress_generators(degree_lhs, lhs_compressed);

  auto rhs_compressed(rhs);
  degree_rhs = compress_generators(degree_rhs, rhs_compressed);

  unsigned degree = degree_lhs * degree_rhs;

  PermSet wreath_product_generators;

  for (unsigned i = 0u; i < degree_rhs; ++i) {
    for (Perm const &perm : lhs_compressed)
      wreath_product_generators.insert(
        perm.shifted(degree_lhs * i).extended(degree));
  }

  for (Perm const &gen_rhs : rhs_compressed) {
    std::vector<std::vector<unsigned>> cycles {gen_rhs.cycles()};
    for (auto &cycle : cycles) {
      for (unsigned &x : cycle)
        x = (x - 1u) * degree_lhs + 1u;
    }

    std::vector<std::vector<unsigned>> shifted_cycles {cycles};

    for (unsigned i = 1u; i < degree_lhs; ++i) {
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
{
  std::vector<int> orbit_indices(degree() + 1u, -1);
  orbit_indices[1] = 0;

  unsigned processed = 1u;

  for (auto i = 1u; i <= degree(); ++i) {
    int orbit_index1 = orbit_indices[i];
    if (orbit_index1 == -1)
      return false;

    for (Perm const &gen : _bsgs.strong_generators()) {
      unsigned j = gen[i];

      int orbit_index2 = orbit_indices[j];
      if (orbit_index2 == -1) {
        orbit_indices[j] = orbit_index1;

        if (++processed == degree())
          return true;
      }
    }
  }

  throw std::logic_error("unreachable");
}

bool PermGroup::contains_element(Perm const &perm) const
{
  assert(perm.degree() == degree() && "element has same degree as group");

  return _bsgs.strips_completely(perm);
}

Perm PermGroup::random_element() const
{
  static std::default_random_engine gen(time(0));

  Perm result(degree());
  for (unsigned i = 0u; i < _bsgs.base_size(); ++i) {
    std::vector<unsigned> orbit = _bsgs.orbit(i);
    std::uniform_int_distribution<> d(0u, orbit.size() - 1u);
    result *= _bsgs.transversal(i, orbit[d(gen)]);
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

      _transversals.emplace_back(transv);
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

std::ostream& operator<<(std::ostream& stream, PermGroup const &pg) {
  stream << pg.bsgs();
  return stream;
}

} // namespace cgtl
