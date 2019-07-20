#include <cassert>
#include <cstddef>
#include <ostream>
#include <set>
#include <vector>

#include "perm.h"
#include "util.h"

/**
 * @file perm.cc
 * @brief Implements `Perm`.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

Perm::Perm(unsigned degree) : _n(degree), _perm(degree)
{
  assert(_n > 0u && "permutation degree > 0");

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

    for (auto i = 1u; i < cycle.size(); ++i) {
      unsigned tmp = cycle[i];
      assert(tmp <= _n && "cycle element <= N");
      _perm[cycle[i - 1u] - 1u] = tmp;
    }

    _perm[cycle.back() - 1u] = cycle[0];

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
  std::vector<unsigned> inverse(_n);

  for (unsigned i = 0u; i < _n; ++i)
    inverse[_perm[i] - 1u] = i + 1u;

  return Perm(inverse);
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
  assert(rhs.degree() == _n && "comparing permutations of equal degree");

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
  assert(rhs.degree() == _n && "multiplying permutations of equal degree");

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

std::vector<std::vector<unsigned>> Perm::cycles() const
{
  std::vector<std::vector<unsigned>> result;

  std::vector<unsigned> cycle;

  std::set<unsigned> done;

  unsigned first, current;
  first = current = 1u;

  while (1) {
    done.insert(current);
    cycle.push_back(current);

    current = (*this)[current];

    if (current == first) {
      if (cycle.size() > 1u)
        result.push_back(cycle);

      cycle.clear();

      if (done.size() == _n)
        return result;

      for (unsigned i = 1u; i <= _n; ++i) {
        if (done.find(i) == done.end()) {
          first = i;
          current = i;
          break;
        }
      }
    }
  }
}

Perm Perm::extended(unsigned degree) const
{
  assert(degree >= _n);

  if (degree == _n)
    return *this;

  std::vector<unsigned> perm(degree);

  for (unsigned i = 0u; i < _n; ++i)
    perm[i] = _perm[i];

  for (unsigned i = _n + 1u; i <= degree; ++i)
    perm[i - 1u] = i;

  return Perm(perm);
}

Perm Perm::normalized(unsigned low, unsigned high) const
{
  std::vector<unsigned> perm_normalized(high - low + 1u);

  for (auto i = low; i <= high; ++i)
    perm_normalized[i - low] = _perm[i - 1u] - low + 1u;

  return Perm(perm_normalized);
}

Perm Perm::shifted(unsigned shift) const
{
  if (shift == 0 || _n == 0u)
    return *this;

  std::vector<unsigned> perm_shifted(_n + shift);

  for (unsigned i = 1u; i <= shift; ++i)
    perm_shifted[i - 1u] = i;

  for (unsigned i = 1u; i <= _n; ++i)
    perm_shifted[i + shift - 1u] = (*this)[i] + shift;

  return Perm(perm_shifted);
}

Perm Perm::restricted(std::vector<unsigned> const &domain) const
{
  // domain must be union of domains of disjoint cycles or behaviour is undefined

  std::vector<unsigned> restricted_perm(_n);
  for (auto i = 1u; i <= _n; ++i)
    restricted_perm[i - 1u] = i;

  for (unsigned x : domain) {
    if (_perm[x - 1u] != x) {
      restricted_perm[x - 1u] = _perm[x - 1u];
    }
  }

  return Perm(restricted_perm);
}

} // namespace cgtl

namespace std
{

std::size_t hash<cgtl::Perm>::operator()(cgtl::Perm const &perm) const
{
  return util::vector_hash(perm._perm);
}

} // namespace std
