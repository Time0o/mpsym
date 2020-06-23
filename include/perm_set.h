#ifndef _GUARD_PERM_SET_H
#define _GUARD_PERM_SET_H

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <ostream>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "dump.h"
#include "perm.h"

namespace mpsym
{

class PermSet
{
  friend std::ostream &operator<<(std::ostream &os, PermSet const &ps);

public:
  typedef std::vector<Perm>::iterator iterator;
  typedef std::vector<Perm>::const_iterator const_iterator;
  typedef std::vector<Perm>::reverse_iterator reverse_iterator;
  typedef std::vector<Perm>::const_reverse_iterator const_reverse_iterator;
  typedef std::vector<Perm>::size_type size_type;

  PermSet()
  {}

  PermSet(std::initializer_list<Perm> perms)
  : PermSet(perms.begin(), perms.end())
  {}

  template<typename IT>
  PermSet(IT b, IT e)
  { insert(b, e); }

  unsigned degree() const
  {
    assert(!empty() && "degree of empty permutation set not defined");
    return _perms[0].degree();
  }

  bool empty() const
  { return _perms.empty(); }

  unsigned size() const
  { return _perms.size(); }

  // TODO: prefer iterators
  PermSet subset(unsigned offs, unsigned sz) const
  {
    assert(offs + sz <= size());
    return PermSet(begin() + offs, begin() + offs + sz);
  }

  Perm const& operator[](unsigned i) const
  {
    assert(i < size());
    return _perms[i];
  }

  Perm & operator[](unsigned i)
  {
    assert(i < size());
    return _perms[i];
  }

  iterator begin() { return _perms.begin(); }
  iterator end() { return _perms.end(); }
  const_iterator begin() const { return _perms.begin(); }
  const_iterator end() const { return _perms.end(); }

  reverse_iterator rbegin() { return _perms.rbegin(); }
  reverse_iterator rend() { return _perms.rend(); }
  const_reverse_iterator rbegin() const { return _perms.rbegin(); }
  const_reverse_iterator rend() const { return _perms.rend(); }

  void insert(Perm const &perm) {
    assert_degree(perm.degree());
    _perms.push_back(perm);
  }

  void insert(Perm &&perm)
  {
    assert_degree(perm.degree());
    _perms.emplace_back(perm); // TODO: forward
  }

  template<typename IT>
  void insert(IT b, IT e)
  {
#ifndef NDEBUG
    for (auto it = b; it != e; ++it)
      assert_degree(it->degree());
#endif
    _perms.insert(_perms.end(), b, e);
  }

  template<typename ...ARGS>
  void emplace(ARGS &&...args)
  {
    _perms.emplace_back(args...); // TODO: forward
    assert_degree(_perms.back().degree());
  }

  size_type erase(Perm const &perm)
  {
    size_type removed = 0u;

    auto it = _perms.begin();
    while (it != _perms.end()) {
      if (*it == perm) {
        it = erase(it);
        ++removed;
      } else {
        ++it;
      }
    }

    return removed;
  }

  template<typename IT>
  IT erase(IT it)
  { return _perms.erase(it); }

  void clear()
  { _perms.clear(); }

  bool trivial() const
  {
    if (empty())
      return true;

    for (auto const &perm : *this) {
      if (!perm.id())
        return false;
    }

    return true;
  }

  bool contains(Perm const &perm) const
  { return std::find(_perms.begin(), _perms.end(), perm) != _perms.end(); }

  unsigned smallest_moved_point() const;
  unsigned largest_moved_point() const;

  void make_unique();
  void make_generating_set();
  void minimize_degree();

  void assert_not_empty() const
  { assert(!empty() && "permutation set not empty"); }

  void assert_degree(unsigned deg) const
  {
#ifndef NDEBUG
    assert((empty() || degree() == deg) && "permutations have correct degree");
#else
    (void)deg;
#endif
  }

  void assert_generating_set() const
  {
    auto unique_perms(unique());

    for (auto const &perm : _perms) {
      assert(unique_perms.find(~perm) != unique_perms.end()
             && "closed under inversion");
    }
  }

private:
  std::unordered_set<Perm> unique() const
  { return std::unordered_set<Perm>(_perms.begin(), _perms.end()); }

  std::vector<Perm> _perms;
};

inline std::ostream &operator<<(std::ostream &os, PermSet const &ps)
{
  std::vector<Perm> perms(ps._perms);
  std::sort(perms.begin(), perms.end());
  os << DUMP_CUSTOM(perms, "{}");
  return os;
}

}

#endif // _GUARD_PERM_SET_H
