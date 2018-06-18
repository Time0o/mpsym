#include <cassert>
#include <utility>
#include <vector>

#include "perm.h"

namespace cgtl
{

PermWord::PermWord(Perm const &perm) : _n(perm.degree())
{
  _perms.push_back(perm);
  _invperms.push_back(~perm);
}

unsigned PermWord::operator[](unsigned const i) const
{
  unsigned result = _perms[0][i];

  for (auto i = 1u; i < _perms.size(); ++i)
    result = _perms[i][result];

  return result;
}

PermWord PermWord::operator~() const
{
  auto n_perms = _perms.size();

  PermWord inverse(*this);
  for (auto i = 0u; i < n_perms; ++i)
    std::swap(inverse._perms[i], inverse._invperms[n_perms - 1u - i]);

  return inverse;
}

PermWord operator*(PermWord const &lhs, PermWord const &rhs)
{
  PermWord result(lhs);
  return result *= rhs;
}

std::ostream& operator<<(std::ostream& stream, PermWord const &pw)
{
  stream << pw.perm();
  return stream;
}

bool PermWord::operator==(PermWord const &rhs) const
{
  for (unsigned i = 1u; i <= degree(); ++i) {
    if ((*this)[i] != rhs[i])
      return false;
  }

  return true;
}

bool PermWord::operator!=(PermWord const &rhs) const
{
  return !((*this) == rhs);
}

PermWord& PermWord::operator*=(PermWord const &rhs)
{
  _perms.insert(_perms.end(), rhs._perms.begin(), rhs._perms.end());
  _invperms.insert(_invperms.end(), rhs._invperms.begin(), rhs._invperms.end());

  return *this;
}

Perm PermWord::perm() const
{
  std::vector<unsigned> perm(_n);
  for (auto i = 0u; i < _n; ++i)
    perm[i] = (*this)[i + 1u];

  return perm;
}

} // namespace cgtl
