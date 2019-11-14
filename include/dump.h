#ifndef _GUARD_DUMP_H
#define _GUARD_DUMP_H

#include <cassert>
#include <cstdint>
#include <initializer_list>
#include <ostream>
#include <set>
#include <type_traits>
#include <unordered_set>
#include <vector>

namespace dump
{

template<typename T>
class Dumper
{
  template<typename U>
  friend std::ostream &operator<<(std::ostream &os, Dumper<U> const &dumper);

  template<typename U>
  struct is_iterable
  {
    typedef uint8_t true_type;
    typedef uint16_t false_type;

    template<typename V> static true_type iterable(
      int,
      typename V::const_iterator = V().begin(),
      typename V::const_iterator = V().end());

    template<typename V> static false_type iterable(...) ;

    enum { value = sizeof(iterable<U>(0)) == sizeof(true_type) };
  };

  template<typename U, template<typename ...> class V>
  struct is_specialization : std::false_type {};

  template<template<typename ...> class V, typename ...W>
  struct is_specialization<V<W...>, V> : std::true_type {};

  template<typename U>
  using is_set = std::integral_constant<
    bool, is_specialization<U, std::set>::value ||
          is_specialization<U, std::unordered_set>::value>;

public:
  Dumper(T const &obj, std::initializer_list<char const *> brackets)
  : _obj(obj),
    _brackets(brackets)
  {}

  void dump(std::ostream &os) const
  { do_dump(os, _obj, 0); }

private:
  template<typename U>
  typename std::enable_if<!is_iterable<U>::value>::type
  do_dump(std::ostream &os, U const &obj, unsigned) const
  { os << obj; }

  template<typename U>
  typename std::enable_if<is_iterable<U>::value>::type
  do_dump(std::ostream &os, U const &obj, unsigned level) const
  {
    os << opening_bracket<U>(level);

    auto it = obj.begin();
    while (it != obj.end()) {
      do_dump(os, *it, level + 1u);

      if (++it != obj.end())
        os << ", ";
    }

    os << closing_bracket<U>(level);
  }

  template<typename U>
  char bracket(bool closing, unsigned level) const
  {
    if (level < _brackets.size())
      return _brackets[level][closing ? 1 : 0];

    return is_set<U>::value ? (closing ? '}' : '{') : (closing ? ']' : '[');
  }

  template<typename U>
  char opening_bracket(unsigned level) const
  { return bracket<U>(false, level); }

  template<typename U>
  char closing_bracket(unsigned level) const
  { return bracket<U>(true, level); }

  T _obj;

  std::vector<char const *> _brackets;
};

template<typename T>
Dumper<T> dump(T const &obj, std::initializer_list<char const *> brackets = {})
{ return Dumper<T>(obj, brackets); }

template<typename T>
std::ostream &operator<<(std::ostream &os, Dumper<T> const &dumper)
{
  dumper.dump(os);
  return os;
}

} // namespace dump

#endif // _GUARD_DUMP_H
