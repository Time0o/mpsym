#ifndef _GUARD_PARTIAL_PERM_H
#define _GUARD_PARTIAL_PERM_H

#include <ostream>
#include <vector>

namespace cgtl
{

class PartialPerm
{
public:
  PartialPerm();
  PartialPerm(std::vector<unsigned> const &pperm);

  unsigned operator[](unsigned const i) const;
  PartialPerm operator~() const;
  bool operator==(PartialPerm const &rhs) const;
  bool operator!=(PartialPerm const &rhs) const;
  PartialPerm& operator*=(PartialPerm const &rhs);

  std::vector<unsigned> image(std::vector<unsigned> const &alpha) const;

  std::vector<unsigned> dom() const { return _dom; }
  unsigned dom_min() const { return _dom_min; }
  unsigned dom_max() const { return _dom_max; }

  std::vector<unsigned> im() const { return _im; }
  unsigned im_min() const { return _im[0]; }
  unsigned im_max() const { return _im[_im.size() - 1u]; }

private:
  std::vector<unsigned> _pperm;
  std::vector<unsigned> _dom, _im;
  unsigned _dom_min, _dom_max;
};

std::ostream& operator<<(std::ostream& stream, PartialPerm const &pperm);
PartialPerm operator*(PartialPerm const &lhs, PartialPerm const &rhs);

} // namespace cgtl

#endif // _GUARD_PARTIAL_PERM_H
