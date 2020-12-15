#ifndef GUARD_PARTIAL_PERM_SET_H
#define GUARD_PARTIAL_PERM_SET_H

#include <cassert>
#include <vector>

#include "partial_perm.hpp"
#include "vector_set.hpp"

namespace mpsym
{

namespace internal
{

class PartialPermSet : public util::VectorSet<PartialPerm, PartialPermSet>
{
public:
  PartialPermSet()
  {}

  PartialPermSet(std::initializer_list<PartialPerm> perms)
  : VectorSet(perms.begin(), perms.end())
  {}

  template<typename IT>
  PartialPermSet(IT first, IT last)
  : VectorSet(first, last)
  {}

  unsigned dom_min() const
  { return get_min(&PartialPerm::dom_min); }

  unsigned dom_max() const
  { return get_max(&PartialPerm::dom_max); }

  unsigned im_min() const
  { return get_min(&PartialPerm::im_min); }

  unsigned im_max() const
  { return get_max(&PartialPerm::im_max); }

private:
  bool check_elem(PartialPerm const &pperm) const override
  { return !pperm.empty(); }

  unsigned get_min(int(PartialPerm::*f)() const) const
  {
    assert(!empty());

    unsigned min = UINT_MAX;
    for (PartialPerm const &pperm : _elems) {
      int tmp = (pperm.*f)();
      assert(tmp != -1);

      min = std::min(min, static_cast<unsigned>(tmp));
    }

    return min;
  }

  unsigned get_max(int(PartialPerm::*f)() const) const
  {
    unsigned max = 0u;
    for (PartialPerm const &pperm : _elems) {
      int tmp = (pperm.*f)();
      assert(tmp != -1);

      max = std::max(max, static_cast<unsigned>(tmp));
    }

    return max;
  }
};

} // namespace internal

} // namespace mpsym

#endif // GUARD_PERM_SET_H
