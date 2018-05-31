#ifndef _GUARD_PERM_H
#define _GUARD_PERM_H

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <vector>

#ifndef NDEBUG
#include <set>
#endif

namespace CGTL
{

class Perm
{
friend Perm operator*(Perm const &lhs, Perm const &rhs) {
  bool left_larger = lhs.size() >= rhs.size();
  Perm result(left_larger ? lhs : rhs);

  if (left_larger) {
    for (unsigned i = 0u; i < rhs.size(); ++i)
      result._perm[i] = lhs[rhs[i + 1u]];
  } else {
    for (unsigned i = 0u; i < rhs.size(); ++i) {
      unsigned tmp = rhs[i + 1u];
      if (tmp <= lhs.size())
        result._perm[i] = lhs[tmp];
    }
  }

  return result;
}

public:
  typedef unsigned size_type;

  Perm(unsigned n = 0u) : _n(n), _perm(n) {
    for (unsigned i = 0u; i < _n; ++i)
       _perm[i] = i + 1u;
  }

  Perm(std::vector<unsigned> perm)
    : _n(*std::max_element(perm.begin(), perm.end())), _perm(perm) {
#ifndef NDEBUG
    assert(("explicit permutation description has correct length",
            perm.size() == _n));

    if (_n == 0u)
      return;

    std::set<unsigned> tmp(perm.begin(), perm.end());

    assert(("explicit permutation description does not contain duplicate elements",
            tmp.size() == perm.size()));

    bool full_range = (*tmp.begin() == 1u) && (*tmp.rbegin() == _n);
    assert(("explicit permutation description contains all elements from 1 to N",
            full_range));
#endif
  }

  Perm(unsigned n, std::vector<std::vector<unsigned>> cycles) : _n(n) {
    assert(_n > 0u);

    if (cycles.size() == 0u) {
      for (unsigned i = 0u; i < _n; ++i)
        _perm.push_back(i + 1u);

    } else if (cycles.size() == 1u) {
      std::vector<unsigned> const &cycle = cycles[0];
#ifndef NDEBUG
      assert(("cycle has plausible length", cycle.size() <= _n));

      std::set<unsigned> tmp(cycle.begin(), cycle.end());
      assert(("cycle does not contain elements > N", *tmp.rbegin() <= _n));
      assert(("cycle does not contain duplicate elements",
              tmp.size() == cycle.size()));
#endif

      for (unsigned i = 0u; i < _n; ++i)
        _perm.push_back(i + 1u);

      for (size_t i = 1u; i < cycle.size(); ++i) {
        unsigned tmp = cycle[i];
        assert(("cycle element <= N", tmp <= _n));
        (*this)[cycle[i - 1u]] = tmp;
      }

      (*this)[cycle.back()] = cycle[0];

    } else {
      Perm result(_n, { cycles.back() });
      for (auto i = cycles.rbegin() + 1; i != cycles.rend(); ++i)
        result = Perm(_n, { *i }) * result;

      _perm = result._perm;
    }
  }

  unsigned const& operator[](unsigned const i) const {
    assert(("permutation index valid", i > 0u && i <= _n));
    return _perm[i - 1u];
  }

  unsigned& operator[](unsigned const i) {
    return const_cast<unsigned&>(static_cast<const Perm*>(this)->operator[](i));
  }

  unsigned size() const { return _n; }

  Perm& extend(unsigned m) {
    assert(("permutation not narrowed", m >= _n));

    for (unsigned i = _n + 1u; i <= m; ++i)
      _perm.push_back(i);

    return *this;
  }

private:
  unsigned _n;
  std::vector<unsigned> _perm;
};

}

#endif // _GUARD_PERM_H
