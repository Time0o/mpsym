#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "dump.h"
#include "orbits.h"
#include "perm.h"
#include "perm_set.h"
#include "pr_randomizer.h"
#include "explicit_transversals.h"
#include "schreier_tree.h"

namespace mpsym
{

void BSGSTransversalsBase::update_schreier_structure(
  unsigned i, unsigned root, unsigned degree, PermSet const &generators)
{
  auto ss(make_schreier_structure(root, degree, generators));

  Orbit::generate(root, generators, ss);

  if (i == _schreier_structures.size())
    _schreier_structures.push_back(ss);
  else
    _schreier_structures[i].swap(ss);
}

void BSGSTransversalsBase::insert_schreier_structure(
  unsigned i, unsigned root, unsigned degree, PermSet const &generators)
{
  _schreier_structures.insert(_schreier_structures.begin() + i, nullptr);

  update_schreier_structure(i, root, degree, generators);
}

BSGS::BSGS(unsigned degree)
: _degree(degree)
{ assert(degree > 0); }

BSGS::BSGS(unsigned degree, PermSet const &generators, Options const *options_)
: _degree(degree)
{
  if (generators.empty() || (generators.size() == 1u && generators[0].id()))
    return;

  assert(degree > 0);

  generators.assert_degree(degree);

  auto options(Options::fill_defaults(options_));

  switch (options.transversals) {
    case Transversals::EXPLICIT:
      _transversals = std::make_shared<BSGSTransversals<ExplicitTransversals>>();
      break;
    case Transversals::SCHREIER_TREES:
      _transversals = std::make_shared<BSGSTransversals<SchreierTree>>();
      break;
    case Transversals::SHALLOW_SCHREIER_TREES:
      throw std::logic_error("TODO");
  }

  DBG(DEBUG) << "=== Constructing BSGS";
  DBG(DEBUG) << "Generators: " << generators;

  if (options.check_altsym && degree > 8u) {
    PrRandomizer pr(generators);

    if (pr.test_symmetric()) {
      construct_symmetric();
    } else if (pr.test_alternating())
      construct_alternating();
    else
      construct_unknown(generators, &options);
  } else {
    construct_unknown(generators, &options);
  }

  DBG(DEBUG) << "==> B = " << _base;
  DBG(DEBUG) << "==> SGS = " << _strong_generators;

  assert(base_size() > 0u);
}

BSGS::order_type BSGS::order() const
{
  order_type res = 1;

  for (unsigned i = 0u; i < base_size(); ++i)
    res *= orbit(i).size();

  return res;
}

PermSet BSGS::strong_generators(unsigned i) const
{
  PermSet ret;
  for (Perm const &sg : strong_generators()) {
    if (sg.stabilizes(_base.begin(), _base.begin() + i))
      ret.insert(sg);
  }

  return ret;
}

Orbit BSGS::orbit(unsigned i) const
{
  auto nodes(schreier_structure(i)->nodes());

  return Orbit(nodes.begin(), nodes.end());
}

Perm BSGS::transversal(unsigned i, unsigned o) const
{ return schreier_structure(i)->transversal(o); }

PermSet BSGS::transversals(unsigned i) const
{
  PermSet transversals;
  for (unsigned o : orbit(i))
    transversals.insert(schreier_structure(i)->transversal(o));

  return transversals;
}

PermSet BSGS::stabilizers(unsigned i) const
{ return schreier_structure(i)->labels(); }

std::pair<Perm, unsigned> BSGS::strip(Perm const &perm, unsigned offs) const
{
  Perm result(perm);

  for (unsigned i = offs; i < base_size(); ++i) {
    unsigned beta = result[base_point(i)];
    if (!schreier_structure(i)->contains(beta))
      return std::make_pair(result, i + 1u);

    result *= ~schreier_structure(i)->transversal(beta);
  }

  return std::make_pair(result, base_size() + 1u);
}

bool BSGS::strips_completely(Perm const &perm) const
{
  auto strip_result(strip(perm));

  return strip_result.first.id() && strip_result.second == base_size() + 1u;
}

void BSGS::extend_base(unsigned bp)
{ _base.push_back(bp); }

void BSGS::extend_base(unsigned bp, unsigned i)
{ _base.insert(_base.begin() + i, bp); }

void BSGS::construct_symmetric()
{
  DBG(DEBUG) << "Group is symmetric";

  if (_degree == 1u)
    return;

  _base.resize(_degree - 1u);
  std::iota(_base.begin(), _base.end(), 1u);

  for (unsigned i = _degree - 1u; i > 0u; --i)
    _strong_generators.insert(Perm(_degree, {{i, _degree}}));

  for (unsigned i = 0u; i < _base.size(); ++i)
    update_schreier_structure(i, _strong_generators.subset(0, _degree - i - 1u));

  _is_symmetric = true;
}

void BSGS::construct_alternating()
{
  DBG(DEBUG) << "Group is alternating";

  if (_degree < 2u)
    return;

  _base.resize(_degree - 2u);
  std::iota(_base.begin(), _base.end(), 1u);

  for (unsigned i = _degree - 2u; i > 0u; --i)
    _strong_generators.insert(Perm(_degree, {{i, _degree - 1u, _degree}}));

  for (unsigned i = 0u; i < _base.size(); ++i)
    update_schreier_structure(i, _strong_generators.subset(0, _degree - i - 2u));

  _is_alternating = true;
}

void BSGS::construct_unknown(PermSet const &generators, Options const *options)
{
  switch (options->construction) {
    case Construction::SCHREIER_SIMS:
      schreier_sims(generators);
      break;
    case Construction::SCHREIER_SIMS_RANDOM:
      schreier_sims_random(generators, options);
      break;
    case Construction::SOLVE:
      solve(generators);
      break;
  }

  if (options->reduce_gens)
    reduce_gens();
}

std::ostream &operator<<(std::ostream &os, BSGS const &bsgs)
{
  os << "BASE: " << DUMP(bsgs._base) << "\n"
     << "SGS: " << DUMP(bsgs._strong_generators);

  return os;
}

} // namespace mpsym
