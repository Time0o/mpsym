#ifndef GUARD_VECTOR_SET_H
#define GUARD_VECTOR_SET_H

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <ostream>
#include <utility>
#include <vector>

#include "dump.hpp"

namespace mpsym
{

namespace util
{

template<typename T, typename T_SET>
class VectorSet
{
  template<typename U, typename U_SET>
  friend std::ostream &operator<<(std::ostream &os, VectorSet<U, U_SET> const &vs);

public:
  using value_type = T;
  using reference = T &;
  using const_reference = T const &;
  using size_type = typename std::vector<T>::size_type;

  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  using reverse_iterator = typename std::vector<T>::reverse_iterator;
  using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

  VectorSet()
  {}

  VectorSet(std::initializer_list<T> elems)
  : VectorSet(elems.begin(), elems.end())
  {}

  template<typename IT>
  VectorSet(IT first, IT last)
  { insert(first, last); }

  bool empty() const
  { return _elems.empty(); }

  size_type size() const
  { return _elems.size(); }

  T_SET subset(size_type offs, size_type sz) const
  {
    assert(offs + sz <= size());
    return T_SET(begin() + offs, begin() + offs + sz);
  }

  const_reference operator[](size_type i) const
  {
    assert(i < size());
    return _elems[i];
  }

  reference operator[](size_type i)
  {
    assert(i < size());
    return _elems[i];
  }

  iterator begin() { return _elems.begin(); }
  iterator end() { return _elems.end(); }
  const_iterator begin() const { return _elems.begin(); }
  const_iterator end() const { return _elems.end(); }

  reverse_iterator rbegin() { return _elems.rbegin(); }
  reverse_iterator rend() { return _elems.rend(); }
  const_reverse_iterator rbegin() const { return _elems.rbegin(); }
  const_reverse_iterator rend() const { return _elems.rend(); }

  void insert(T const &elem)
  {
    assert_elem(elem);
    _elems.push_back(elem);
  }

  void insert(T &&elem)
  {
    _elems.emplace_back(std::forward<T>(elem));
    assert_elem(_elems.back());
  }

  template<typename IT>
  void insert(IT first, IT last)
  { _elems.insert(_elems.end(), first, last); }

  void resize(size_type n)
  { _elems.resize(n); }

  void resize(size_type n, value_type const &value)
  { _elems.resize(n, value); }

  template<typename ...ARGS>
  void emplace(ARGS &&...args)
  {
    _elems.emplace_back(std::forward<ARGS>(args)...);
    assert_elem(_elems.back());
  }

  size_type erase(T const &elem)
  {
    size_type removed = 0u;

    auto it = _elems.begin();
    while (it != _elems.end()) {
      if (*it == elem) {
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
  { return _elems.erase(it); }

  void clear()
  { _elems.clear(); }

  bool contains(T const &elem) const
  { return std::find(_elems.begin(), _elems.end(), elem) != _elems.end(); }

protected:
  std::vector<T> _elems;

private:
  void assert_elem(T const &elem) const
  { assert(check_elem(elem)); }

  virtual bool check_elem(T const &elem) const
  { return true; }
};

template<typename T, typename T_SET>
inline std::ostream &operator<<(std::ostream &os, VectorSet<T, T_SET> const &vs)
{
  std::vector<T> elems(vs._elems);
  std::sort(elems.begin(), elems.end());
  os << DUMP_CUSTOM(elems, "{}");
  return os;
}

} // namespace util

} // namespace mpsym

#endif // GUARD_VECTOR_SET_H
