#ifndef _GUARD_PERM_GROUP_H
#define _GUARD_PERM_GROUP_H

#include <map>
#include <vector>

#include "bsgs.h"
#include "perm.h"

namespace cgtl
{

class PermGroup
{
public:
  class const_iterator
  {
  public:
    const_iterator() : _end(true) {};
    const_iterator(BSGS const &bsgs);

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
    bool _end;

    std::vector<std::vector<Perm>> _transversals;
    std::vector<Perm> _current_factors;
    Perm _current_result;
  };

  PermGroup(unsigned degree, std::vector<Perm> const &generators);

  static PermGroup symmetric(unsigned degree);
  static PermGroup cyclic(unsigned degree);
  static PermGroup alternating(unsigned degree);

  const_iterator begin() const { return const_iterator(_bsgs); }
  const_iterator end() const { return const_iterator(); }

  unsigned degree() const { return _n; }
  unsigned order() const { return _order; }
  BSGS bsgs() const { return _bsgs; }

  bool is_element(Perm const &perm);
  Perm random_element();

private:
  unsigned _n;
  unsigned _order;
  BSGS _bsgs;
};

} // namespace cgtl

#endif // _GUARD_PERM_GROUP_H
