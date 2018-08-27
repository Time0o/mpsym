#ifndef _GUARD_PERM_H
#define _GUARD_PERM_H

#include <cstddef>
#include <functional>
#include <ostream>
#include <vector>

#include "boost/container_hash/hash.hpp"

namespace cgtl
{

class Perm;
class PermWord;

} // namespace cgtl

namespace std
{

template<>
struct hash<cgtl::Perm>
{
  std::size_t operator()(cgtl::Perm const &perm) const;
};

template<>
struct hash<cgtl::PermWord>
{
  std::size_t operator()(cgtl::PermWord const &perm_word) const;
};

} // namespace std

namespace cgtl
{

class Perm
{
friend std::size_t std::hash<Perm>::operator()(Perm const &perm) const;

public:
  Perm(unsigned degree = 1);
  Perm(std::vector<unsigned> const &perm);
  Perm(unsigned n, std::vector<std::vector<unsigned>> const &cycles);

  unsigned const& operator[](unsigned const i) const;
  Perm operator~() const;
  bool operator==(Perm const &rhs) const;
  bool operator!=(Perm const &rhs) const;
  Perm& operator*=(Perm const &rhs);

  unsigned degree() const { return _n; }
  bool id() const;

  Perm shifted(unsigned low, unsigned high) const;

  Perm restricted(
    std::vector<unsigned> const &domain, bool *id = nullptr) const;

private:
  unsigned _n;
  std::vector<unsigned> _perm;
};

class PermWord
{
friend std::size_t std::hash<PermWord>::operator()(PermWord const &permWord) const;

public:
  PermWord(Perm const &perm);
  PermWord(unsigned degree = 1) : PermWord(Perm(degree)) {};
  PermWord(std::vector<unsigned> const &perm) : PermWord(Perm(perm)) {};
  PermWord(unsigned n, std::vector<std::vector<unsigned>> const &cycles)
    : PermWord(Perm(n, cycles)) {};

  unsigned operator[](unsigned const i) const;
  PermWord operator~() const;
  bool operator==(PermWord const &rhs) const;
  bool operator!=(PermWord const &rhs) const;
  PermWord& operator*=(PermWord const &rhs);

  unsigned degree() const { return _n; }
  bool id() const;
  Perm perm() const;

private:
  unsigned _n;
  std::vector<Perm> _perms;
  std::vector<Perm> _invperms;
};

std::ostream& operator<<(std::ostream& stream, Perm const &perm);
Perm operator*(Perm const &lhs, Perm const &rhs);

std::ostream& operator<<(std::ostream& stream, PermWord const &pw);
PermWord operator*(PermWord const &lhs, PermWord const &rhs);

inline std::size_t hash_value(Perm const &perm) { return std::hash<Perm>()(perm); }

} // namespace cgtl

#endif // _GUARD_PERM_H
