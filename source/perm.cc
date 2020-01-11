#include <algorithm>
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

Perm::Perm(unsigned deg)
: _degree(deg),
  _perm(deg + 1u)
{
  assert(degree() > 0u && "permutation degree > 0");

  for (unsigned i = 1u; i <= degree(); ++i)
    _perm[i] = i;
}

Perm::Perm(std::vector<unsigned> const &perm)
: _degree(*std::max_element(perm.begin(), perm.end())),
  _perm(perm.size() + 1u)
{
  assert(perm.size() == degree() &&
    "explicit permutation description has correct length");

#ifndef NDEBUG
  if (degree() == 0u)
    return;

  std::set<unsigned> tmp(perm.begin(), perm.end());
  bool full_range = (*tmp.begin() == 1u) && (*tmp.rbegin() == degree());
#endif

  assert(tmp.size() == degree() &&
    "explicit permutation description does not contain duplicate elements");

  assert(full_range &&
    "explicit permutation description contains all elements from 1 to N");

  std::copy(perm.begin(), perm.end(), _perm.begin() + 1u);
}

Perm::Perm(unsigned n, std::vector<std::vector<unsigned>> const &cycles)
: Perm(n)
{
  assert(degree() > 0u);

  if (cycles.size() == 0u)
    return;

  if (cycles.size() == 1u) {
    std::vector<unsigned> const &cycle = cycles[0];

    assert(cycle.size() <= degree() && "cycle has plausible length");

#ifndef NDEBUG
    std::set<unsigned> tmp(cycle.begin(), cycle.end());
#endif

    assert(*tmp.rbegin() <= degree() &&
      "cycle does not contain elements > N");
    assert(tmp.size() == cycle.size() &&
      "cycle does not contain duplicate elements");

    for (auto i = 1u; i < cycle.size(); ++i) {
      unsigned tmp = cycle[i];
      assert(tmp <= degree() && "cycle element <= N");
      _perm[cycle[i - 1u]] = tmp;
    }

    _perm[cycle.back()] = cycle[0];

  } else {
    for (auto i = cycles.begin(); i != cycles.end(); ++i)
      (*this) *= Perm(degree(), {*i});
  }
}

unsigned const& Perm::operator[](unsigned const i) const
{
  assert(i > 0u && i <= degree() && "permutation index valid");
  return _perm[i];
}

Perm Perm::operator~() const
{
  std::vector<unsigned> inverse(degree());

  for (unsigned i = 1u; i <= degree(); ++i)
    inverse[(*this)[i] - 1u] = i;

  return Perm(inverse);
}

Perm operator*(Perm const &lhs, Perm const &rhs)
{
  Perm result(lhs);
  return result *= rhs;
}

std::ostream &operator<<(std::ostream &os, const Perm &perm)
{
  std::set<unsigned> done;

  unsigned first, current;
  first = current = 1u;

  bool id = true;

  std::vector<unsigned> cycle;

  for (;;) {
    done.insert(current);
    cycle.push_back(current);

    current = perm[current];

    if (current == first) {
      if (cycle.size() > 1u) {
        id = false;

        os << "(" << cycle[0];
        for (auto i = 1u; i < cycle.size(); ++i)
          os << " " << cycle[i];
        os << ")";

      }

      cycle.clear();

      if (done.size() == perm.degree()) {
        if (id)
          os << "()";
        return os;
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
  assert(rhs.degree() == degree() && "comparing permutations of equal degree");

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
  assert(rhs.degree() == degree() && "multiplying permutations of equal degree");

  for (unsigned i = 1u; i <= rhs.degree(); ++i)
    _perm[i] = rhs[(*this)[i]];

  return *this;
}

bool Perm::id() const
{
  for (unsigned i = 1u; i <= degree(); ++i) {
    if ((*this)[i] != i)
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

  for (;;) {
    done.insert(current);
    cycle.push_back(current);

    current = (*this)[current];

    if (current == first) {
      if (cycle.size() > 1u)
        result.push_back(cycle);

      cycle.clear();

      if (done.size() == degree())
        return result;

      for (unsigned i = 1u; i <= degree(); ++i) {
        if (done.find(i) == done.end()) {
          first = i;
          current = i;
          break;
        }
      }
    }
  }
}

Perm Perm::extended(unsigned deg) const
{
  assert(deg >= degree());

  if (deg == degree())
    return *this;

  std::vector<unsigned> perm(deg);

  for (unsigned i = 0u; i < degree(); ++i)
    perm[i] = (*this)[i + 1u];

  for (unsigned i = degree() + 1u; i <= deg; ++i)
    perm[i - 1u] = i;

  return Perm(perm);
}

Perm Perm::normalized(unsigned low, unsigned high) const
{
  std::vector<unsigned> perm_normalized(high - low + 1u);

  for (auto i = low; i <= high; ++i)
    perm_normalized[i - low] = (*this)[i] - low + 1u;

  return Perm(perm_normalized);
}

Perm Perm::shifted(unsigned shift) const
{
  if (shift == 0 || degree() == 0u)
    return *this;

  std::vector<unsigned> perm_shifted(degree() + shift);

  for (unsigned i = 1u; i <= shift; ++i)
    perm_shifted[i - 1u] = i;

  for (unsigned i = 1u; i <= degree(); ++i)
    perm_shifted[i + shift - 1u] = (*this)[i] + shift;

  return Perm(perm_shifted);
}

} // namespace cgtl

namespace std
{

std::size_t hash<cgtl::Perm>::operator()(cgtl::Perm const &perm) const
{ return util::container_hash(perm._perm.begin() + 1u, perm._perm.end()); }

} // namespace std
