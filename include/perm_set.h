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
  typedef std::vector<Perm>::const_iterator const_iterator;

  PermSet()
  {}

  PermSet(std::initializer_list<Perm> perms)
  : PermSet(perms.begin(), perms.end())
  {}

  template<typename IT>
  PermSet(IT b, IT e)
  {
    insert(b, e);
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

  const_iterator begin() const { return _perms.begin(); }
  const_iterator end() const { return _perms.end(); }

  void insert(Perm const &perm) {
    assert_degree(perm.degree());
    _perms.push_back(perm);
  }

  void insert(Perm &&perm) {
    assert_degree(perm.degree());
    _perms.emplace_back(perm); // TODO: forward
  }

  template<typename IT>
  void insert(IT b, IT e) {
#ifndef NDEBUG
    for (auto it = b; it != e; ++it)
      assert_degree(it->degree());
#endif
    _perms.insert(_perms.end(), b, e);
  }

  template<typename ...ARGS>
  void emplace(ARGS &&...args) {
    _perms.emplace_back(args...); // TODO: forward
    assert_degree(_perms.back().degree());
  }

  template<typename IT>
  IT erase(IT it) { return _perms.erase(it); }

  void clear() { _perms.clear(); }

  void make_unique();

  void minimize_degree();

  void assert_not_empty() const {
    assert(!empty() && "permutation set not empty");
  }

  void assert_degree(unsigned deg) const {
#ifndef NDEBUG
    assert((empty() || degree() == deg) && "permutations have correct degree");
#else
    (void)deg;
#endif
  }

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
