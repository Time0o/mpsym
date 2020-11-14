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
  PartialPerm(std::vector<int> const &dom,
              std::vector<int> const &im);
  PartialPerm(std::vector<int> const &pperm);

  int operator[](int i) const;
  PartialPerm operator~() const;
  bool operator==(PartialPerm const &rhs) const;
  bool operator!=(PartialPerm const &rhs) const;
  PartialPerm& operator*=(PartialPerm const &rhs);

  static PartialPerm from_perm(Perm const &perm);
  Perm to_perm(unsigned degree) const;

  std::vector<int> dom() const
  { return _dom; }

  int dom_min() const
  { return _dom.empty() ? -1 : _dom[0]; }

  int dom_max() const
  { return _dom.empty() ? -1 : _dom.back(); }

  std::vector<int> im() const
  { return _im; }

  int im_min() const
  { return _im.empty() ? -1 : _im[0]; }

  int im_max() const
  { return _im.empty() ? -1 : _im.back(); }

  bool empty() const
  { return _pperm.empty(); }

  bool id() const
  { return _id; }

  template<typename IT>
  PartialPerm restricted(IT first, IT last) const
  {
    if (first == last)
      return PartialPerm();

    std::vector<int> pperm_restricted(dom_max(), -1);

    for (IT it = first; it != last; ++it) {
      int x = *it;

      if (x < dom_min() || x > dom_max())
        continue;

      unsigned y = (*this)[x];

      if (y != 0u)
        pperm_restricted[x] = y;
    }

    if (!pperm_restricted.empty()) {
      int dom_max_restricted = dom_max();

      while (pperm_restricted[dom_max_restricted - 1] == -1)
        --dom_max_restricted;

      pperm_restricted.resize(dom_max_restricted);
    }

    return PartialPerm(pperm_restricted);
  }

  template<template<typename ...> class T, typename IT>
  T<unsigned> image(IT first, IT last) const
  {
    std::set<int> pperm_image;

    for (IT it = first; it != last; ++it) {
      int x = *it;

      if (x < dom_min() || x > dom_max())
        continue;

      int y = _pperm[x];

      if (y != -1)
        pperm_image.insert(y);
    }

    return T<int>(pperm_image.begin(), pperm_image.end());
  }

private:
  std::vector<int> _pperm, _dom, _im;
  bool _id;
};

std::ostream &operator<<(std::ostream &os, PartialPerm const &pperm);

PartialPerm operator*(PartialPerm const &lhs, PartialPerm const &rhs);

} // namespace internal

} // namespace mpsym

#endif // GUARD_PARTIAL_PERM_H
