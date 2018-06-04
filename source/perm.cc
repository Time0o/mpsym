#include <cassert>
#include <set>
#include <ostream>
#include <vector>

#include "perm.h"

namespace cgtl
{

Perm::Perm(unsigned degree) : _n(degree), _perm(degree)
{
  assert(_n > 0u && "Permutation degree > 0");

  for (unsigned i = 0u; i < _n; ++i)
     _perm[i] = i + 1u;
}

Perm::Perm(std::vector<unsigned> const &perm) : _perm(perm)
{
  _n = 0u;
  for (unsigned u : perm) {
    if (u > _n)
      _n = u;
  }

  assert(perm.size() == _n &&
    "explicit permutation description has correct length");

#ifndef NDEBUG
  if (_n == 0u)
    return;

  std::set<unsigned> tmp(perm.begin(), perm.end());
  bool full_range = (*tmp.begin() == 1u) && (*tmp.rbegin() == _n);
#endif

  assert(tmp.size() == perm.size() &&
    "explicit permutation description does not contain duplicate elements");

  assert(full_range &&
    "explicit permutation description contains all elements from 1 to N");
}

Perm::Perm(unsigned n, std::vector<std::vector<unsigned>> const &cycles)
 : Perm(n)
{
  assert(_n > 0u);

  if (cycles.size() == 0u)
    return;

  if (cycles.size() == 1u) {
    std::vector<unsigned> const &cycle = cycles[0];

    assert(cycle.size() <= _n && "cycle has plausible length");

#ifndef NDEBUG
    std::set<unsigned> tmp(cycle.begin(), cycle.end());
#endif

    assert(*tmp.rbegin() <= _n &&
      "cycle does not contain elements > N");
    assert(tmp.size() == cycle.size() &&
      "cycle does not contain duplicate elements");

    for (size_t i = 1u; i < cycle.size(); ++i) {
      unsigned tmp = cycle[i];
      assert(tmp <= _n && "cycle element <= N");
      (*this)._perm[cycle[i - 1u] - 1u] = tmp;
    }

    (*this)._perm[cycle.back() - 1u] = cycle[0];

  } else {
    for (auto i = cycles.begin(); i != cycles.end(); ++i)
      (*this) *= Perm(_n, {*i});
  }
}

unsigned const& Perm::operator[](unsigned const i) const
{
  assert(i > 0u && i <= _n && "permutation index valid");
  return _perm[i - 1u];
}

Perm Perm::operator~() const
{
  std::vector<unsigned> perm_inverted(_perm);

  for (unsigned i = 0u; i < _n; ++i)
    perm_inverted[_perm[i] - 1u] = i + 1u;

  return Perm(perm_inverted);
}

Perm operator*(Perm const &lhs, Perm const &rhs)
{
  Perm result(lhs);
  return result *= rhs;
}

std::ostream& operator<<(std::ostream& stream, const Perm &perm)
{
  std::set<unsigned> done;

  unsigned first, current;
  first = current = 1u;

  bool id = true;

  std::vector<unsigned> cycle;

  while (1) {
    done.insert(current);
    cycle.push_back(current);

    current = perm[current];

    if (current == first) {
      if (cycle.size() > 1u) {
        id = false;

        stream << "(" << cycle[0];
        for (auto i = 1u; i < cycle.size(); ++i)
          stream << " " << cycle[i];
        stream << ")";

      }

      cycle.clear();

      if (done.size() == perm.degree()) {
        if (id)
          stream << "()";
        return stream;
      }

      for (unsigned i = 1u; i <= perm.degree(); ++i) {
        if (done.find(i) == done.end()) {
          first = i;
          current = i;
          break;
        }
      }
    }
  }
}

bool Perm::operator==(Perm const &rhs) const
{
  for (unsigned i = 1u; i <= degree(); ++i) {
    if ((*this)[i] != rhs[i])
      return false;
  }

  return true;
}

bool Perm::operator!=(Perm const &rhs) const
{
  return !((*this) == rhs);
}

Perm& Perm::operator*=(Perm const &rhs)
{
  for (unsigned i = 0u; i < rhs.degree(); ++i)
    _perm[i] = rhs[(*this)[i + 1u]];

  return *this;
}

bool Perm::id() const
{
  for (unsigned i = 0u; i < _n; ++i) {
    if (_perm[i] != i + 1u)
      return false;
  }
  return true;
}

} // namespace cgtl
