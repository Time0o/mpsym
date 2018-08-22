#ifndef _GUARD_PERM_GROUP_H
#define _GUARD_PERM_GROUP_H

#include <map>
#include <vector>

#include "bsgs.h"
#include "perm.h"
#include "schreier_sims.h"

namespace cgtl
{

class PermGroup
{
public:
  class const_iterator
  {
  public:
    const_iterator() : _end(true) {};
    const_iterator(PermGroup const &pg);

    const_iterator operator++();
    const_iterator operator++(int) { next_state(); return *this; }
    Perm const & operator*() const { return _current_result; }
    Perm const * operator->() const { return &_current_result; }
    bool operator==(const_iterator const &rhs) const;
    bool operator!=(const_iterator const &rhs) const { return !((*this) == rhs); }

  private:
    void next_state();
    void update_result();

    std::vector<unsigned> _state;
    bool _trivial;
    bool _end;

    std::vector<std::vector<Perm>> _transversals;
    std::vector<Perm> _current_factors;
    Perm _current_result;
  };

  PermGroup() : PermGroup(1, {}) {}
  PermGroup(unsigned degree, std::vector<Perm> const &generators,
    SchreierSims::Variant schreier_sims_method = SchreierSims::SIMPLE);

  static PermGroup symmetric(unsigned degree);
  static PermGroup cyclic(unsigned degree);
  static PermGroup alternating(unsigned degree);

  const_iterator begin() const { return const_iterator(*this); }
  const_iterator end() const { return const_iterator(); }

  unsigned degree() const { return _n; }
  unsigned order() const { return _order; }
  BSGS bsgs() const { return _bsgs; }
  bool trivial() const { return _bsgs.trivial(); }
  bool transitive() const;

  bool is_element(Perm const &perm) const;
  Perm random_element() const;

  std::vector<PermGroup> disjoint_decomposition(
    bool complete = true, bool disjoint_orbit_optimization = false) const;

private:
  std::vector<PermGroup> disjoint_decomposition_complete(
    bool disjoint_orbit_optimization = true) const;

  std::vector<PermGroup> disjoint_decomposition_incomplete() const;

  unsigned _n;
  unsigned _order;
  BSGS _bsgs;
};

} // namespace cgtl

#endif // _GUARD_PERM_GROUP_H
