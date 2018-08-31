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
  if (pperm.empty())
    return;

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
  if (pperm.dom().empty()) {
    stream << "()";
    return stream;
  }

  std::vector<std::vector<unsigned>> chains;
  std::vector<std::vector<unsigned>> cycles;

  unsigned first, current;
  first = current = pperm.dom_min();

  std::set<unsigned> done;
  std::vector<unsigned> current_chain;

  for (;;) {
    done.insert(current);
    current_chain.push_back(current);

    current = pperm[current];

    bool end_of_chain = current == 0u || current > pperm.dom_max();
    bool end_of_cycle = current == first;

    if (!(end_of_chain || end_of_cycle))
      continue;

    if (end_of_chain) {
      if (current != 0u)
        current_chain.push_back(current);

      if (current_chain.size() > 1u)
        chains.push_back(current_chain);

    } else if (end_of_cycle)
      cycles.push_back(current_chain);

    current_chain.clear();

    bool next = false;
    for (unsigned i = pperm.dom_min(); i <= pperm.dom_max(); ++i) {
      if (done.find(i) == done.end()) {
        next = true;
        first = i;
        current = i;
        break;
      }
    }

    if (!next)
      break;
  }

  auto sort_pred_size = [](std::vector<unsigned> const &a,
                           std::vector<unsigned> const &b) {
    return a.size() > b.size();
  };

  std::sort(chains.begin(), chains.end(), sort_pred_size);

  std::vector<std::vector<unsigned>> unique_chains;

  for (auto const &chain : chains) {
    bool superfluous = false;
    for (auto const &other : unique_chains) {
      if (other.size() < chain.size())
        continue;
      auto it = std::search(other.begin(), other.end(),
                            chain.begin(), chain.end());

      if (it != other.end()) {
        superfluous = true;
        break;
      }
    }

    if (!superfluous)
      unique_chains.push_back(chain);
  }

  auto sort_pred_init_elem = [](std::vector<unsigned> const &a,
                                std::vector<unsigned> const &b) {
    return a[0] < b[0];
  };

  std::sort(unique_chains.begin(), unique_chains.end(), sort_pred_init_elem);

  for (auto const &chain : unique_chains) {
     stream << '[' << chain[0];
     for (auto i = 1u; i < chain.size(); ++i)
       stream << ' ' << chain[i];
     stream << ']';
  }

  for (auto const &cycle : cycles) {
    stream << '(' << cycle[0];
    for (auto i = 1u; i < cycle.size(); ++i)
      stream << ' ' << cycle[i];
    stream << ')';
  }

  return stream;
}

} // namespace cgtl
