#include <utility>
#include <vector>

#include "bsgs.h"
#include "perm.h"
#include "schreier_tree.h"
#include "schreier_sims.h"

namespace cgtl
{

BSGS::BSGS(std::vector<unsigned> const &base, std::vector<Perm> const &generators)
  : _base(base), _sgs(generators)
{
  SchreierSims::schreier_sims(_base, _sgs, _schreier_trees);

  for (auto const &st : _schreier_trees)
    _base_elems.push_back(BaseElem(st));
}

std::pair<Perm, unsigned> BSGS::strip(Perm const &perm) const
{
  return SchreierSims::strip(perm, _base, _schreier_trees);
}

std::vector<Perm> BSGS::BaseElem::transversals() const
{
  std::vector<Perm> transversals;
  for (unsigned o : orbit())
    transversals.push_back(transversal(o));

  return transversals;
}

} // namespace cgtl
