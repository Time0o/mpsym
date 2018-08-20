#include <algorithm>
#include <cassert>
#include <ostream>
#include <set>
#include <vector>

#include "partial_perm.h"

namespace cgtl
{

PartialPerm::PartialPerm(std::vector<unsigned> const &pperm) : _pperm(pperm)
{
  for (auto i = 1u; i <= _pperm.size(); ++i) {
    unsigned im = _pperm[i - 1u];

    if (im != 0u) {
      _dom.push_back(i);
      _im.push_back(im);
    }
  }

  _dom_min = _dom[0];
  _dom_max = _dom[_dom.size() - 1u];

  std::sort(_im.begin(), _im.end());

  assert(std::unique(_im.begin(), _im.end()) == _im.end() &&
    "partial permutation image does not contain duplicates");
}

unsigned PartialPerm::operator[](unsigned const i) const
{
  if (i < _dom_min || i > _dom_max)
    return 0u;

  return _pperm[i - 1u];
}

std::vector<unsigned> PartialPerm::image(
  std::vector<unsigned> const &alpha) const
{
  std::vector<unsigned> res;

  for (unsigned x : alpha) {
    if (x < _dom_min || x > _dom_max)
      continue;

    unsigned y = _pperm[x - 1u];
    if (y != 0u)
      res.push_back(y);
  }

  return res;
}

std::ostream& operator<<(std::ostream& stream, PartialPerm const &pperm)
{
  std::set<unsigned> done;

  unsigned first, current;
  first = current = pperm.dom_min();

  std::vector<unsigned> chain;

  for (;;) {
    done.insert(current);
    chain.push_back(current);

    current = pperm[current];

    if (current == 0u || current == first) {

      if (current == 0u && chain.size() != 0u) {
        stream << '[' << chain[0];
        for (auto i = 1u; i < chain.size(); ++i)
          stream << ' ' << chain[i];
        stream << ']';

      } else if (current == first) {
        stream << '(' << chain[0];
        for (auto i = 1u; i < chain.size(); ++i)
          stream << ' ' << chain[i];
        stream << ')';
      }

      chain.clear();

      if (done.size() == pperm.dom_max() - pperm.dom_min() + 1u)
        return stream;

      for (unsigned i = pperm.dom_min(); i <= pperm.dom_max(); ++i) {
        if (done.find(i) == done.end()) {
          first = i;
          current = i;
          break;
        }
      }
    }
  }
}

} // namespace cgtl
