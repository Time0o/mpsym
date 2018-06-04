#ifndef _GUARD_PERM_H
#define _GUARD_PERM_H

#include <vector>

namespace cgtl
{

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

std::ostream& operator<<(std::ostream& stream, const Perm &perm);
Perm operator*(Perm const &lhs, Perm const &rhs);

} // namespace cgtl

#endif // _GUARD_PERM_H
