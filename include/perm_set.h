#ifndef _GUARD_PERM_SET_H
#define _GUARD_PERM_SET_H

#include <cassert>
#include <initializer_list>
#include <ostream>
#include <vector>

#include "perm.h"

namespace cgtl
{

class PermSet
{
public:
  PermSet()
  {}

  PermSet(std::initializer_list<Perm> const &perms)
  : PermSet(perms.begin(), perms.end())
  {}

  // TODO: remove
  PermSet(std::vector<Perm> const &perms)
  : PermSet(perms.begin(), perms.end())
  {}

  template<typename IT>
  PermSet(IT b, IT e)
  : _perms(b, e)
  {
#ifndef NDEBUG
    for (auto i = 1u; i < _perms.size(); ++i)
      assert(_perms[i].degree() == degree());
#endif
  }

  unsigned degree() const {
    assert(!empty() && "degree of empty permutation set not defined");
    return _perms[0].degree();
  }

  bool empty() const { return _perms.empty(); }
  unsigned size() const { return _perms.size(); }

  Perm const& operator[](unsigned i) const {
    assert(i < size());
    return _perms[i];
  };

  Perm & operator[](unsigned i) {
    assert(i < size());
    return _perms[i];
  }

  std::vector<Perm>::const_iterator begin() const { return _perms.begin(); }
  std::vector<Perm>::const_iterator end() const { return _perms.end(); }

  void add(Perm const &perm) { _perms.push_back(perm); }
  void add(Perm &&perm) { _perms.emplace_back(perm); }

  template<typename ...ARGS>
  void add(ARGS &&...args) {  _perms.emplace_back(args...); }

  void assert_not_empty() const {
    assert(!empty() && "permutation set not empty");
  }

  void assert_degree(unsigned n) const {
    assert((empty() || degree() == n) && "permutations have correct degree");
  }

  // TODO: remove
  std::vector<Perm> vect() const { return _perms; }

private:
  std::vector<Perm> _perms;
};

inline std::ostream& operator<<(std::ostream& stream, PermSet const &perm_set)
{
  stream << "[";
  if (!perm_set.empty()) {
    stream << perm_set[0];
    for (auto i = 1u; i < perm_set.size(); ++i)
      stream << ", " << perm_set[i];
  }
  stream << "]";

  return stream;
}

}

#endif // _GUARD_PERM_SET_H
