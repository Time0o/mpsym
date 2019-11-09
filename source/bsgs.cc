#include <algorithm>
#include <cassert>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "orbits.h"
#include "perm.h"
#include "perm_set.h"
#include "explicit_transversals.h"
#include "schreier_tree.h"

namespace cgtl
{

BSGS::BSGS(unsigned degree,
           PermSet const &generators,
           Construction construction,
           Transversals transversals)
: _degree(degree),
  _transversals(transversals)
{
  generators.assert_degree(degree);

  if (generators.empty() || (generators.size() == 1u && generators[0].id()))
    return;

  switch (construction) {
    case Construction::SCHREIER_SIMS:
      schreier_sims(generators);
      break;
    case Construction::SCHREIER_SIMS_RANDOM:
      schreier_sims_random(generators);
      break;
    case Construction::SOLVE:
      solve(generators);
      break;
    case Construction::AUTO:
      schreier_sims(generators); // TODO
      break;
  }

  assert(base_size() > 0u);

  reduce_gens();
}

std::vector<unsigned> BSGS::orbit(unsigned i) const
{
  return _schreier_structures[i]->nodes();
}

Perm BSGS::transversal(unsigned i, unsigned o) const
{
  return _schreier_structures[i]->transversal(o);
}

PermSet BSGS::transversals(unsigned i) const
{
  PermSet transversals;
  for (unsigned o : orbit(i))
    transversals.insert(_schreier_structures[i]->transversal(o));

  return transversals;
}

PermSet BSGS::stabilizers(unsigned i) const
{
  return _schreier_structures[i]->labels();
}

std::pair<Perm, unsigned> BSGS::strip(Perm const &perm, unsigned offs) const
{
  Perm result(perm);

  for (unsigned i = offs; i < base_size(); ++i) {
    unsigned beta = result[base_point(i)];
    if (!_schreier_structures[i]->contains(beta))
      return std::make_pair(result, i + 1u);

    result *= ~_schreier_structures[i]->transversal(beta);
  }

  return std::make_pair(result, base_size() + 1u);
}

bool BSGS::strips_completely(Perm const &perm) const
{
  auto strip_result(strip(perm));

  return strip_result.first.id() && strip_result.second == base_size() + 1u;
}

void BSGS::extend_base(unsigned bp)
{
  _base.push_back(bp);
}

void BSGS::update_schreier_structure(unsigned i, PermSet const &generators)
{
  std::shared_ptr<SchreierStructure> ss;

  auto d = degree();
  auto r = base_point(i);
  auto gens(generators);

  switch (_transversals) {
    case Transversals::EXPLICIT:
      ss = std::make_shared<ExplicitTransversals>(d, r, gens);
      break;
    case Transversals::SCHREIER_TREES:
    case Transversals::AUTO:
      ss = std::make_shared<SchreierTree>(d, r, gens);
      break;
    case Transversals::SHALLOW_SCHREIER_TREES:
      throw std::logic_error("TODO");
  }

  orbit_of(r, gens, ss.get());

  if (_schreier_structures.size() < base_size()) {
    assert(i == _schreier_structures.size());
    _schreier_structures.emplace_back(std::move(ss));
  } else {
    _schreier_structures[i].swap(ss);
  }
}

std::ostream& operator<<(std::ostream& stream, BSGS const &bsgs)
{
  stream << "BASE: [";

  for (auto i = 0u; i < bsgs.base_size(); ++i) {
    stream << bsgs.base_point(i);
    if (i < bsgs.base_size() - 1u)
      stream << ", ";
  }

  stream << "]; SGS: [";

  for (auto i = 0u; i < bsgs._strong_generators.size(); ++i) {
    stream << bsgs._strong_generators[i];
    if (i < bsgs._strong_generators.size() - 1u)
      stream << ", ";
  }

  stream << ']';
  return stream;
}

} // namespace cgtl
