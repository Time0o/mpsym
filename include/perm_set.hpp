#ifndef GUARD_PERM_SET_H
#define GUARD_PERM_SET_H

#include <cassert>
#include <vector>

#include "perm.hpp"
#include "vector_set.hpp"

namespace mpsym
{

namespace internal
{

class PermSet : public util::VectorSet<Perm, PermSet>
{
public:
  PermSet()
  {}

  PermSet(std::initializer_list<Perm> perms)
  : VectorSet(perms.begin(), perms.end())
  {}

  template<typename IT>
  PermSet(IT first, IT last)
  : VectorSet(first, last)
  {}

  unsigned degree() const
  { return _elems[0].degree(); }

  bool trivial() const;

  unsigned smallest_moved_point() const;
  unsigned largest_moved_point() const;
  std::vector<unsigned> support() const;

  bool has_inverses() const;
  PermSet with_inverses() const;
  void insert_inverses();

  void make_unique();
  void minimize_degree();

  void assert_degree(unsigned deg) const
  {
#ifndef NDEBUG
    assert((empty() || degree() == deg) && "permutations have correct degree");
#else
    (void)deg;
#endif
  }

  void assert_inverses() const
  { assert(has_inverses() && "closed under inversion"); }

private:
  bool check_elem(Perm const &perm) const override
  { return empty() || perm.degree() == degree(); }
};

} // namespace internal

} // namespace mpsym

#endif // GUARD_PERM_SET_H
