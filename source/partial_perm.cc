#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>
#include <ostream>
#include <set>
#include <vector>

#include "partial_perm.h"
#include "perm.h"
#include "util.h"

/**
 * @file partial_perm.cc
 * @author Timo Nicolai
 *
 * @brief Implements `PartialPerm`.
 */

namespace mpsym
{

namespace internal
{

PartialPerm::PartialPerm(unsigned degree)
: _pperm(degree),
  _id(true)
{
  std::iota(_pperm.begin(), _pperm.end(), 1u);

  _dom = _pperm;
  _im = _pperm;
}

PartialPerm::PartialPerm(std::vector<unsigned> const &dom,
                         std::vector<unsigned> const &im)
: _dom(dom),
  _im(im),
  _id(true)
{
  assert(dom.size() == im.size() &&
         "partial permutation domain and image have same dimension");

  if (dom.empty())
    return;

  _pperm = std::vector<unsigned>(*std::max_element(_dom.begin(), _dom.end()), 0u);

  for (auto i = 0u; i < _dom.size(); ++i) {
    unsigned x = _dom[i];
    unsigned y = _im[i];

    assert(x != 0u && y != 0u
           && "partial permutation domain and image do not contain 0 elements");

    if (_id && y != x)
      _id = false;

    _pperm[x - 1u] = y;
  }

  std::sort(_dom.begin(), _dom.end());
  std::sort(_im.begin(), _im.end());

  assert(std::unique(_dom.begin(), _dom.end()) == _dom.end() &&
         "partial permutation domain does not contain duplicate elements");

  assert(std::unique(_im.begin(), _im.end()) == _im.end() &&
         "partial permutation image does not contain duplicate elements");
}

PartialPerm::PartialPerm(std::vector<unsigned> const &pperm)
: _pperm(pperm),
  _id(true)
{
  if (pperm.empty())
    return;

  for (auto i = 1u; i <= _pperm.size(); ++i) {
    unsigned im = _pperm[i - 1u];

    if (_id && im != 0u && im != i)
      _id = false;

    if (im != 0u) {
      _dom.push_back(i);
      _im.push_back(im);
    }
  }

  std::sort(_im.begin(), _im.end());

  assert(std::unique(_im.begin(), _im.end()) == _im.end() &&
    "partial permutation image does not contain duplicates");
}

unsigned PartialPerm::operator[](unsigned const i) const
{
  assert(i > 0u && i <= _pperm.size() && "partial permutation index valid");
  return _pperm[i - 1u];
}

PartialPerm PartialPerm::operator~() const
{
  PartialPerm res;

  std::vector<unsigned> inverse(im_max(), 0u);
  for (unsigned x : _dom) {
    unsigned y = (*this)[x];
    inverse[y - 1u] = x;
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
     os << '[' << chain[0];
     for (auto i = 1u; i < chain.size(); ++i)
       os << ' ' << chain[i];
     os << ']';
  }

  for (auto const &cycle : cycles) {
    os << '(' << cycle[0];
    for (auto i = 1u; i < cycle.size(); ++i)
      os << ' ' << cycle[i];
    os << ')';
  }

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
    unsigned y = (*this)[x];

    unsigned z;
    if (y < rhs.dom_min() || y > rhs.dom_max())
      z = 0u;
    else
      z = rhs[y];

    if (z != 0u) {
      dom_new.push_back(x);
      im_new.push_back(z);
    }

    _pperm[x - 1u] = z;

    if (_id && x != z)
      _id = false;
  }

  _dom = dom_new;
  _im = im_new;

  std::sort(_im.begin(), _im.end());

  decltype(_pperm.size()) reduce = 0u;
  for (auto i = _pperm.size() - 1u; i > 0u; --i) {
    if (_pperm[i] != 0u)
      break;

    ++reduce;
  }

  _pperm.resize(_pperm.size() - reduce);

  return *this;
}

PartialPerm PartialPerm::from_perm(Perm const &perm)
{
  std::vector<unsigned> pperm(perm.degree());

  for (unsigned i = 1u; i <= perm.degree(); ++i)
    pperm[i - 1u] = perm[i];

  return PartialPerm(pperm);
}

Perm PartialPerm::to_perm(unsigned degree) const
{
  std::vector<unsigned> perm(degree);

  unsigned i;
  for (i = 1u; i <= degree; ++i) {
    if (i < dom_min()) {
      perm[i - 1u] = i;
    } else if (i > dom_max()) {
      break;
    } else {
      unsigned im = _pperm[i - 1u];
      if (im == 0u)
        perm[i - 1u] = i;
      else if (i > 0u) {
#ifndef NDEBUG
        auto it = std::find(perm.begin(), perm.begin() + i, im);
        assert(it == perm.begin() + i &&
               "partial permutation does not contain chain within domain");
#endif
        perm[i - 1u] = im;
      }
    }
  }

  while (i <= degree) {
    perm[i - 1u] = i;
    ++i;
  }

  return Perm(perm);
}

} // namespace internal

} // namespace mpsym

namespace std
{

std::size_t hash<mpsym::PartialPerm>::operator()(
  mpsym::PartialPerm const &pperm) const
{ return util::container_hash(pperm._pperm.begin(), pperm._pperm.end()); }

} // namespace std
