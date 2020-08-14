#ifndef GUARD_PARTIAL_PERM_H
#define GUARD_PARTIAL_PERM_H

#include <cstddef>
#include <ostream>
#include <set>
#include <vector>

namespace mpsym
{

namespace internal
{

class PartialPerm;
class Perm;

} // namespace internal

} // namespace mpsym

namespace std
{

template<>
struct hash<mpsym::internal::PartialPerm>
{
  std::size_t operator()(mpsym::internal::PartialPerm const &pperm) const;
};

} // namespace std

namespace mpsym
{

namespace internal
{

class PartialPerm
{
friend std::size_t std::hash<PartialPerm>::operator()(
  PartialPerm const &perm) const;

public:
  PartialPerm(unsigned degree = 0u);
  PartialPerm(std::vector<unsigned> const &dom,
              std::vector<unsigned> const &im);
  PartialPerm(std::vector<unsigned> const &pperm);

  unsigned operator[](unsigned const i) const;
  PartialPerm operator~() const;
  bool operator==(PartialPerm const &rhs) const;
  bool operator!=(PartialPerm const &rhs) const;
  PartialPerm& operator*=(PartialPerm const &rhs);

  static PartialPerm from_perm(Perm const &perm);
  Perm to_perm(unsigned degree) const;

  std::vector<unsigned> dom() const
  { return _dom; }

  unsigned dom_min() const
  { return _dom.empty() ? 0u : _dom[0]; }

  unsigned dom_max() const
  { return _dom.empty() ? 0u : _dom.back(); }

  std::vector<unsigned> im() const
  { return _im; }

  unsigned im_min() const
  { return _im.empty() ? 0u : _im[0]; }

  unsigned im_max() const
  { return _im.empty() ? 0u : _im.back(); }

  bool empty() const
  { return _pperm.empty(); }

  bool id() const
  { return _id; }

  template<typename IT>
  PartialPerm restricted(IT first, IT last) const
  {
    if (first == last)
      return PartialPerm();

    std::vector<unsigned> pperm_restricted(dom_max(), 0u);

    for (IT it = first; it != last; ++it) {
      unsigned x = *it;

      if (x < dom_min() || x > dom_max())
        continue;

      unsigned y = (*this)[x];

      if (y != 0u)
        pperm_restricted[x - 1u] = y;
    }

    if (!pperm_restricted.empty()) {
      unsigned dom_max_restricted = dom_max();

      while (pperm_restricted[dom_max_restricted - 1u] == 0u)
        --dom_max_restricted;

      pperm_restricted.resize(dom_max_restricted);
    }

    return PartialPerm(pperm_restricted);
  }

  template<template<typename ...> class T, typename IT>
  T<unsigned> image(IT first, IT last) const
  {
    std::set<unsigned> res;

    for (IT it = first; it != last; ++it) {
      unsigned x = *it;

      if (x < dom_min() || x > dom_max())
        continue;

      unsigned y = _pperm[x - 1u];

      if (y != 0u)
        res.insert(y);
    }

    return T<unsigned>(res.begin(), res.end());
  }

private:
  std::vector<unsigned> _pperm, _dom, _im;
  bool _id;
};

std::ostream &operator<<(std::ostream &os, PartialPerm const &pperm);

PartialPerm operator*(PartialPerm const &lhs, PartialPerm const &rhs);

} // namespace internal

} // namespace mpsym

#endif // GUARD_PARTIAL_PERM_H
