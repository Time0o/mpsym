#include <map>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "perm.h"
#include "schreier_sims.h"

namespace cgtl
{

BSGS::BSGS(std::vector<Perm> const &generators,
  schreier_sims::Variant schreier_sims_method) : _sgs(generators)
{
  if (_sgs.size() > 0u) {

    switch (schreier_sims_method) {
      case schreier_sims::RANDOM:
        schreier_sims::schreier_sims_random(_base, _sgs, _schreier_trees);
        break;
      default:
        schreier_sims::schreier_sims(_base, _sgs, _schreier_trees);
    }

    for (auto const &st : _schreier_trees)
      _base_elems.push_back(BaseElem(st));
  }
}

std::vector<Perm> BSGS::stabilizers(unsigned i) const
{
  std::vector<Perm> res;
  for (Perm const &perm : _sgs) {
    bool stabilizes = true;

    for (auto j = 1u; j <= i; ++j) {
      unsigned b = _base[j - 1u];
      if (perm[b] != b) {
        stabilizes = false;
        break;
      }
    }

    if (stabilizes)
      res.push_back(perm);
  }

  return res;
}

std::vector<std::vector<unsigned>> BSGS::orbits() const
{
  std::vector<std::vector<unsigned>> result;

  for (auto const &st : _schreier_trees)
    result.push_back(st.nodes());

  return result;
}

std::pair<Perm, unsigned> BSGS::strip(Perm const &perm) const
{
  return schreier_sims::strip(perm, _base, _schreier_trees);
}

std::vector<Perm> BSGS::BaseElem::transversals() const
{
  std::vector<Perm> transversals;
  for (unsigned o : orbit())
    transversals.push_back(transversal(o));

  return transversals;
}

} // namespace cgtl
