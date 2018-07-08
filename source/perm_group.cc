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
    _current_result *= _current_factors[j];
}

} // namespace cgtl
