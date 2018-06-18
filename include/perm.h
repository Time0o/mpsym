#ifndef _GUARD_PERM_H
#define _GUARD_PERM_H

#include <ostream>
#include <vector>

namespace cgtl
{

class Perm;
class PermWord;

class Perm
{
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

private:
  unsigned _n;
  std::vector<unsigned> _perm;
};

class PermWord
{
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

} // namespace cgtl

#endif // _GUARD_PERM_H
