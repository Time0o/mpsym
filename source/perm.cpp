#include <algorithm>
#include <cassert>
#include <iterator>
#include <numeric>
#include <ostream>
#include <set>
#include <vector>

#include "dump.hpp"
#include "perm.hpp"
#include "util.hpp"

namespace mpsym
{

namespace internal
{

Perm::Perm(unsigned deg)
: _degree(deg),
  _perm(deg)
{
  assert(degree() > 0u);

  std::iota(_perm.begin(), _perm.end(), 0u);
}

Perm::Perm(std::vector<unsigned> const &perm)
{
  _perm = perm;

  assert(!_perm.empty());

  _degree = *std::max_element(_perm.begin(), _perm.end()) + 1u;

  assert(_perm.size() == degree());

#ifndef NDEBUG
  std::set<unsigned> domain(_perm.begin(), _perm.end());

  assert(domain.size() == degree());
  assert(*domain.begin() == 0u);
  assert(*domain.rbegin() == degree() - 1u);
#endif
}

Perm::Perm(unsigned deg, std::vector<std::vector<unsigned>> const &cycles)
: Perm(deg)
{
  if (cycles.size() == 0u)
    return;

  if (cycles.size() == 1u) {
    std::vector<unsigned> const &cycle = cycles[0];

    assert(cycle.size() <= degree());

#ifndef NDEBUG
    std::set<unsigned> cycle_domain(cycle.begin(), cycle.end());

    assert(*cycle_domain.rbegin() < degree());
    assert(cycle_domain.size() == cycle.size());
#endif

    for (auto i = 1u; i < cycle.size(); ++i) {
      assert(cycle[i] < degree());
      _perm[cycle[i - 1u]] = cycle[i];
    }

    _perm[cycle.back()] = cycle[0];

  } else {
    for (auto i = cycles.begin(); i != cycles.end(); ++i)
      (*this) *= Perm(degree(), {*i});
  }
}

unsigned const& Perm::operator[](unsigned i) const
{
  assert(i < degree());
  return _perm[i];
}

Perm Perm::operator~() const
{
  std::vector<unsigned> inverse(degree());

  for (unsigned i = 0u; i < degree(); ++i)
    inverse[(*this)[i]] = i;

  return Perm(inverse);
}

std::ostream &operator<<(std::ostream &os, const Perm &perm)
{
  if (perm.id()) {
    os << "()";
  } else {
    for (auto const &cycle : perm.cycles())
      os << DUMP_CUSTOM(cycle, "()");
  }

  return os;
}

bool Perm::operator==(Perm const &rhs) const
{
  assert(rhs.degree() == degree());

  for (unsigned i = 0u; i < degree(); ++i) {
    if ((*this)[i] != rhs[i])
      return false;
  }

  return true;
}

bool Perm::operator<(Perm const &rhs) const
{ return cycles() < rhs.cycles(); }

Perm& Perm::operator*=(Perm const &rhs)
{
  assert(rhs.degree() == degree());

  for (unsigned i = 0u; i < rhs.degree(); ++i)
    _perm[i] = rhs[(*this)[i]];

  return *this;
}

bool Perm::id() const
{
  for (unsigned i = 0u; i < degree(); ++i) {
    if ((*this)[i] != i)
      return false;
  }
  return true;
}

bool Perm::even() const
{
  bool odd = false;

  for (auto const &cycle : cycles()) {
    if (cycle.size() % 2 == 0)
      odd = !odd;
  }

  return !odd;
}

std::vector<std::vector<unsigned>> Perm::cycles() const
{
  std::vector<std::vector<unsigned>> result;

  std::vector<unsigned> cycle;

  std::set<unsigned> done;

  unsigned first, current;
  first = current = 0u;

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

      for (unsigned i = 0u; i < degree(); ++i) {
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
    perm[i] = (*this)[i];

  for (unsigned i = degree(); i < deg; ++i)
    perm[i] = i;

  return Perm(perm);
}

Perm Perm::normalized(unsigned low, unsigned high) const
{
  std::vector<unsigned> perm_normalized(high - low + 1u);

  for (auto i = low; i <= high; ++i)
    perm_normalized[i - low] = (*this)[i] - low;

  return Perm(perm_normalized);
}

Perm Perm::shifted(unsigned shift) const
{
  if (shift == 0)
    return *this;

  std::vector<unsigned> perm_shifted(degree() + shift);

  for (unsigned i = 0u; i < shift; ++i)
    perm_shifted[i] = i;

  for (unsigned i = 0u; i < degree(); ++i)
    perm_shifted[i + shift] = (*this)[i] + shift;

  return Perm(perm_shifted);
}

} // namespace internal

} // namespace mpsym

namespace std
{

std::size_t hash<mpsym::internal::Perm>::operator()(
  mpsym::internal::Perm const &perm) const
{ return mpsym::util::container_hash(perm._perm.begin() + 1u, perm._perm.end()); }

} // namespace std
