#include <algorithm>
#include <cassert>
#include <numeric>
#include <ostream>
#include <set>
#include <vector>

#include "dump.hpp"
#include "partial_perm.hpp"
#include "perm.hpp"
#include "util.hpp"

namespace mpsym
{

namespace internal
{

PartialPerm::PartialPerm(unsigned degree)
: _pperm(degree),
  _dom(degree),
  _im(degree),
  _id(true)
{
  std::iota(_pperm.begin(), _pperm.end(), 0);
  std::iota(_dom.begin(), _dom.end(), 0u);
  std::iota(_im.begin(), _im.end(), 0u);
}

PartialPerm::PartialPerm(std::vector<unsigned> const &dom,
                         std::vector<unsigned> const &im)
: _dom(dom),
  _im(im),
  _id(true)
{
  assert(_dom.size() == _im.size() &&
         "partial permutation domain and image have same dimension");

  if (_dom.empty())
    return;

  unsigned degree = *std::max_element(_dom.begin(), _dom.end()) + 1;

  assert(degree > 0);

  _pperm = std::vector<int>(degree, -1);

  for (auto i = 0u; i < _dom.size(); ++i) {
    unsigned x = _dom[i];
    unsigned y = _im[i];

    if (_id && y != x)
      _id = false;

    _pperm[x] = y;
  }

  std::sort(_dom.begin(), _dom.end());
  std::sort(_im.begin(), _im.end());

  assert(std::unique(_dom.begin(), _dom.end()) == _dom.end() &&
         "partial permutation domain does not contain duplicate elements");

  assert(std::unique(_im.begin(), _im.end()) == _im.end() &&
         "partial permutation image does not contain duplicate elements");
}

PartialPerm::PartialPerm(std::vector<int> const &pperm)
: _pperm(pperm),
  _id(true)
{
  if (pperm.empty())
    return;

  for (int x = 0; x < static_cast<int>(_pperm.size()); ++x) {
    int y = _pperm[x];

    assert(y >= -1);

    if (y != -1) {
      _dom.push_back(x);
      _im.push_back(y);

      if (_id && y != x)
        _id = false;
    }
  }

  std::sort(_im.begin(), _im.end());

  assert(std::unique(_im.begin(), _im.end()) == _im.end() &&
    "partial permutation image does not contain duplicates");
}

int PartialPerm::operator[](int i) const
{
  assert(i < _pperm.size());
  return _pperm[i];
}

PartialPerm PartialPerm::operator~() const
{
  PartialPerm res;

  std::vector<int> inverse(im_max() + 1, -1);
  for (int x : _dom) {
    int y = (*this)[x];
    inverse[y] = x;
  }

  res._pperm = inverse;
  res._dom = _im;
  res._im = _dom;

  return res;
}

PartialPerm operator*(PartialPerm const &lhs, PartialPerm const &rhs)
{
  PartialPerm result(lhs);
  return result *= rhs;
}

std::ostream &operator<<(std::ostream &os, PartialPerm const &pperm)
{
  if (pperm.dom().empty()) {
    os << "()";
    return os;
  }

  std::vector<std::vector<int>> chains;
  std::vector<std::vector<int>> cycles;

  int first, current;
  first = current = pperm.dom_min();

  std::set<int> done;
  std::vector<int> current_chain;

  for (;;) {
    done.insert(current);
    current_chain.push_back(current);

    current = pperm[current];

    bool end_of_chain = current == -1 || current > pperm.dom_max();
    bool end_of_cycle = current == first;

    if (!end_of_chain && !end_of_cycle)
      continue;

    if (end_of_chain) {
      if (current != -1)
        current_chain.push_back(current);

      if (current_chain.size() > 1u)
        chains.push_back(current_chain);

    } else if (end_of_cycle)
      cycles.push_back(current_chain);

    current_chain.clear();

    bool next = false;
    for (int x = pperm.dom_min(); x <= pperm.dom_max(); ++x) {
      if (done.find(x) != done.end())
        continue;

      next = true;

      first = current = x;

      break;
    }

    if (!next)
      break;
  }

  auto sort_pred_size = [](std::vector<int> const &a,
                           std::vector<int> const &b) {
    return a.size() > b.size();
  };

  std::sort(chains.begin(), chains.end(), sort_pred_size);

  std::vector<std::vector<int>> unique_chains;

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

  auto sort_pred_init_elem = [](std::vector<int> const &a,
                                std::vector<int> const &b) {
    return a[0] < b[0];
  };

  std::sort(unique_chains.begin(), unique_chains.end(), sort_pred_init_elem);

  for (auto const &chain : unique_chains)
    os << DUMP_CUSTOM(chain, "[]");

  for (auto const &cycle : cycles)
    os << DUMP_CUSTOM(cycle, "()");

  return os;
}

bool PartialPerm::operator==(PartialPerm const &rhs) const
{
  if (rhs.dom_min() != dom_min() || rhs.dom_max() != dom_max())
    return false;

  return _pperm == rhs._pperm;
}

bool PartialPerm::operator!=(PartialPerm const &rhs) const
{ return !(*this == rhs); }

PartialPerm& PartialPerm::operator*=(PartialPerm const &rhs)
{
  if (_dom.empty())
    return *this;

  if (rhs.dom().empty()) {
    _pperm.clear();
    _dom.clear();
    _im.clear();

    _id = true;

    return *this;
  }

  std::vector<unsigned> dom_new;
  std::vector<unsigned> im_new;

  _id = true;

  for (unsigned x : _dom) {
    int y = (*this)[static_cast<int>(x)];

    int z;
    if (y == -1 || y < rhs.dom_min() || y > rhs.dom_max()) {
      z = -1;
    } else {
      z = rhs[y];
    }

    if (z != -1) {
      dom_new.push_back(static_cast<unsigned>(x));
      im_new.push_back(static_cast<unsigned>(z));

      if (_id && z != x)
        _id = false;
    }

    _pperm[x] = z;
  }

  _dom = dom_new;
  _im = im_new;

  std::sort(_im.begin(), _im.end());

  decltype(_pperm.size()) reduce = 0u;
  for (auto i = _pperm.size() - 1u; i > 0u; --i) {
    if (_pperm[i] != -1)
      break;

    ++reduce;
  }

  _pperm.resize(_pperm.size() - reduce);

  return *this;
}

PartialPerm PartialPerm::from_perm(Perm const &perm)
{
  std::vector<int> pperm(perm.degree());

  for (unsigned i = 0u; i < perm.degree(); ++i)
    pperm[i] = static_cast<int>(perm[i]);

  return PartialPerm(pperm);
}

Perm PartialPerm::to_perm(unsigned degree) const
{
  if (_dom.empty())
    return Perm(degree);

  std::vector<unsigned> perm(degree);

  unsigned x;
  for (x = 0u; x < degree; ++x) {
    if (x < static_cast<unsigned>(dom_min())) {
      perm[x] = x;
    } else if (x > static_cast<unsigned>(dom_max())) {
      break;
    } else {
      int y = _pperm[x];
      if (y == -1)
        perm[x] = x;
      else
        perm[x] = y;
    }
  }

  while (x < degree) {
    perm[x] = x;
    ++x;
  }

  return Perm(perm);
}

int PartialPerm::dom_min() const
{ return _dom.empty() ? -1 : static_cast<int>(_dom[0]); }

int PartialPerm::dom_max() const
{ return _dom.empty() ? -1 : static_cast<int>(_dom.back()); }

int PartialPerm::im_min() const
{ return _im.empty() ? -1 : static_cast<int>(_im[0]); }

int PartialPerm::im_max() const
{ return _im.empty() ? -1 : static_cast<int>(_im.back()); }

} // namespace internal

} // namespace mpsym

namespace std
{

std::size_t hash<mpsym::internal::PartialPerm>::operator()(
  mpsym::internal::PartialPerm const &pperm) const
{ return mpsym::util::container_hash(pperm._pperm.begin(), pperm._pperm.end()); }

} // namespace std
