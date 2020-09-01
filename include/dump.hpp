#ifndef GUARD_DUMP_H
#define GUARD_DUMP_H

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <ostream>
#include <set>
#include <type_traits>
#include <unordered_set>
#include <vector>

namespace mpsym
{

namespace internal
{

namespace dump
{

template<typename U>
class is_dumpable
{
  template<typename V>
  static auto dumpable(V const *v) -> decltype(std::cout << *v);

  static auto dumpable(...) -> std::false_type;

public:
  enum { value = !std::is_same<decltype(dumpable((U*)0)), std::false_type>::value };
};

template<typename U>
class is_iterable
{
  template<typename V>
  static auto iterable(V const *v) -> decltype(V().begin());

  static auto iterable(...) -> std::false_type;

public:
  enum { value = !std::is_same<decltype(iterable((U*)0)), std::false_type>::value };
};

template<typename U, template<typename ...> class V>
struct is_specialization : std::false_type {};

template<template<typename ...> class V, typename ...W>
struct is_specialization<V<W...>, V> : std::true_type {};

template<typename U>
using is_set = std::integral_constant<
  bool, is_specialization<U, std::set>::value ||
        is_specialization<U, std::unordered_set>::value>;

template<typename T>
class Dumper
{
  template<typename U>
  friend std::ostream &operator<<(std::ostream &os, Dumper<U> const &dumper);

public:
  Dumper(T const &obj, std::initializer_list<char const *> brackets)
  : _obj(obj),
    _brackets(brackets)
  {}

  void dump(std::ostream &os) const
  { do_dump(os, _obj, 0); }

private:
  template<typename U>
  typename std::enable_if<is_dumpable<U>::value || !is_iterable<U>::value>::type
  do_dump(std::ostream &os, U const &obj, unsigned) const
  { os << obj; }

  template<typename U>
  typename std::enable_if<!is_dumpable<U>::value && is_iterable<U>::value>::type
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
dump::Dumper<T> make_dumper(T const &obj,
                            std::initializer_list<char const *> brackets = {})
{ return dump::Dumper<T>(obj, brackets); }

template<typename T, typename FUNC>
using transformed_type =
  typename std::result_of<FUNC(decltype(*std::declval<T>().begin()))>::type;

template<typename T, typename FUNC>
using transformed_container = typename std::conditional<
                                is_set<T>::value,
                                std::set<transformed_type<T, FUNC>>,
                                std::vector<transformed_type<T, FUNC>>
                              >::type;

template<typename T, typename FUNC>
dump::Dumper<transformed_container<T, FUNC>>
transform_and_make_dumper(
  T const &obj,
  FUNC &&trans,
  std::initializer_list<char const *> brackets = {})
{
  static_assert(is_iterable<T>::value, "is iterable");

  transformed_container<T, FUNC> cont;
  std::transform(obj.begin(), obj.end(), std::back_inserter(cont), trans);
  return dump::Dumper<decltype(cont)>(cont, brackets);
}

template<typename T>
std::ostream &operator<<(std::ostream &os, Dumper<T> const &dumper)
{
  dumper.dump(os);
  return os;
}

} // namespace dump

} // namespace internal

} // namespace mpsym

#define DUMP_NS ::mpsym::internal::dump

#define DUMP(obj) DUMP_NS :: make_dumper(obj)
#define TRANSFORM_AND_DUMP(obj, func) DUMP_NS :: transform_and_make_dumper(obj, func)

#define DUMP_CUSTOM(obj, ...) DUMP_NS :: make_dumper(obj, { __VA_ARGS__ })
#define TRANFORM_AND_DUMP_CUSTOM(obj, func, ...) DUMP_NS :: transform_and_make_dumper(obj, func, { __VA_ARGS__ })

#endif // GUARD_DUMP_H
