#include <algorithm>
#include <cassert>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "dump.h"
#include "orbits.h"
#include "perm.h"
#include "perm_set.h"
#include "explicit_transversals.h"
#include "schreier_tree.h"

namespace cgtl
{

void BSGSTransversalsBase::update_schreier_structure(
  unsigned i, unsigned root, PermSet const &generators)
{
  auto ss(make_schreier_structure(root, generators));

  Orbit::generate(root, generators, ss);

  if (i == _schreier_structures.size())
    _schreier_structures.push_back(ss);
  else
    _schreier_structures[i].swap(ss);
}

void BSGSTransversalsBase::insert_schreier_structure(
  unsigned i, unsigned root, PermSet const &generators)
{
  _schreier_structures.insert(_schreier_structures.begin() + i, nullptr);

  update_schreier_structure(i, root, generators);
}

BSGS::BSGS(unsigned degree)
: _degree(degree)
{}

BSGS::BSGS(unsigned degree,
           PermSet const &generators,
           Construction construction,
           Transversals transversals)
: _degree(degree)
{
  generators.assert_degree(degree);

  switch (transversals) {
    case Transversals::AUTO:
    case Transversals::EXPLICIT:
      _transversals = std::make_shared<BSGSTransversals<ExplicitTransversals>>();
      break;
    case Transversals::SCHREIER_TREES:
      _transversals = std::make_shared<BSGSTransversals<SchreierTree>>();
      break;
    case Transversals::SHALLOW_SCHREIER_TREES:
      throw std::logic_error("TODO");
  }

  switch (construction) {
    case Construction::AUTO:
    case Construction::SCHREIER_SIMS:
      schreier_sims(generators);
      break;
    case Construction::SCHREIER_SIMS_RANDOM:
      schreier_sims_random(generators);
      break;
    case Construction::SOLVE:
      solve(generators);
      break;
  }

  assert(base_size() > 0u);

  reduce_gens();
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

std::ostream &operator<<(std::ostream &os, BSGS const &bsgs)
{
  os << "BASE: " << DUMP(bsgs._base)
     << "; SGS: " << DUMP(bsgs._strong_generators);

  return os;
}

} // namespace cgtl
