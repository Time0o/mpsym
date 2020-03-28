#ifndef _GUARD_PERM_WORD_H
#define _GUARD_PERM_WORD_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "perm.h"
#include "perm_set.h"

/**
 * @file perm_word.h
 * @brief Defines `PermWord`.
 *
 * @author Timo Nicolai
 */

namespace mpsym
{

class PermWord;

inline std::size_t hash_value(Perm const &perm)
{
  return std::hash<mpsym::Perm>()(perm);
}

} // namespace mpsym

namespace std
{

template<>
struct hash<mpsym::PermWord>
{
  std::size_t operator()(mpsym::PermWord const &perm_word) const;
};

} // namespace std

namespace mpsym
{

class PermWord
{
friend std::size_t std::hash<PermWord>::operator()(PermWord const &permWord) const;

public:
  PermWord(Perm const &perm);

  PermWord(unsigned degree = 1)
  : PermWord(Perm(degree))
  {};

  PermWord(std::vector<unsigned> const &perm)
  : PermWord(Perm(perm))
  {};

  PermWord(unsigned degree, std::vector<std::vector<unsigned>> const &cycles)
  : PermWord(Perm(degree, cycles))
  {};

  unsigned operator[](unsigned const i) const;
  PermWord operator~() const;
  bool operator==(PermWord const &rhs) const;
  bool operator!=(PermWord const &rhs) const;
  PermWord& operator*=(PermWord const &rhs);

  unsigned degree() const { return _degree; }
  bool id() const;
  Perm perm() const;

private:
  unsigned _degree;
  PermSet _perms;
  PermSet _invperms;
};

std::ostream &operator<<(std::ostream &os, PermWord const &pw);
PermWord operator*(PermWord const &lhs, PermWord const &rhs);

} // namespace mpsym

#endif // _GUARD_PERM_WORD_H
