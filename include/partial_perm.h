#ifndef _GUARD_PARTIAL_PERM_H
#define _GUARD_PARTIAL_PERM_H

#include <ostream>
#include <vector>

#include "perm.h"

namespace cgtl
{

class PartialPerm
{
public:
  PartialPerm(unsigned degree = 0);
  PartialPerm(std::vector<unsigned> const &dom,
              std::vector<unsigned> const &im);
  PartialPerm(std::vector<unsigned> const &pperm);

  static PartialPerm id(std::vector<unsigned> const &dom);
  static PartialPerm from_perm(Perm const &perm);

  unsigned operator[](unsigned const i) const;
  PartialPerm operator~() const;
  bool operator==(PartialPerm const &rhs) const;
  bool operator!=(PartialPerm const &rhs) const;
  PartialPerm& operator*=(PartialPerm const &rhs);

  std::vector<unsigned> dom() const { return _dom; }
  unsigned dom_min() const { return _dom_min; }
  unsigned dom_max() const { return _dom_max; }

  std::vector<unsigned> im() const { return _im; }
  unsigned im_min() const { return _im_min; }
  unsigned im_max() const { return _im_max; }

  bool empty() const { return _pperm.empty(); }
  bool id() const { return _id; }

  PartialPerm restricted(std::vector<unsigned> const &domain) const;
  Perm to_perm(unsigned degree) const;

  std::vector<unsigned> image(std::vector<unsigned> const &alpha) const;

private:
  void update_limits();

  std::vector<unsigned> _pperm;
  std::vector<unsigned> _dom, _im;
  unsigned _dom_min, _dom_max;
  unsigned _im_min, _im_max;
  bool _id;
};

std::ostream& operator<<(std::ostream& stream, PartialPerm const &pperm);
PartialPerm operator*(PartialPerm const &lhs, PartialPerm const &rhs);

} // namespace cgtl

#endif // _GUARD_PARTIAL_PERM_H
