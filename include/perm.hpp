#ifndef GUARD_PERM_H
#define GUARD_PERM_H

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <vector>

#include <boost/operators.hpp>

namespace mpsym
{

namespace internal
{

class Perm;

} // namespace internal

} // namespace mpsym

namespace std
{

template<>
struct hash<mpsym::internal::Perm>
{
  std::size_t operator()(mpsym::internal::Perm const &perm) const;
};

} // namespace std

namespace mpsym
{

namespace internal
{

class Perm : boost::operators<Perm>
{
friend std::size_t std::hash<Perm>::operator()(Perm const &perm) const;

public:
  explicit Perm(unsigned degree = 1);

  explicit Perm(std::vector<unsigned> const &perm);

  Perm(unsigned degree, std::vector<std::vector<unsigned>> const &cycles);

  unsigned const& operator[](unsigned const x) const;
  Perm operator~() const;
  bool operator==(Perm const &rhs) const;
  bool operator<(Perm const &rhs) const;
  Perm& operator*=(Perm const &rhs);

  unsigned degree() const { return _degree; }
  bool id() const;
  bool even() const;

  template<typename IT>
  bool stabilizes(IT first, IT last) const
  {
    for (IT it = first; it != last; ++it) {
       unsigned x = *it;

       if ((*this)[x] != x)
         return false;
    }

    return true;
  }

  Perm extended(unsigned degree) const;
  Perm normalized(unsigned low, unsigned high) const;
  Perm shifted(unsigned shift) const;

  template<typename IT>
  Perm restricted(IT first, IT last) const
  {
    std::vector<std::vector<unsigned>> restricted_cycles;

    for (auto const &cycle : cycles()) {
      bool restrict = false;

      for (unsigned x : cycle) {
        if (std::find(first, last, x) == last)
          restrict = true;
      }

      if (!restrict)
        restricted_cycles.push_back(cycle);
    }

    return Perm(degree(), restricted_cycles);
  }

  std::vector<unsigned> vect() const { return _perm; }
  std::vector<std::vector<unsigned>> cycles() const;

private:
  unsigned _degree;
  std::vector<unsigned> _perm;
};

std::ostream &operator<<(std::ostream &os, Perm const &perm);

} // namespace internal

} // namespace mpsym

#endif // GUARD_PERM_H
