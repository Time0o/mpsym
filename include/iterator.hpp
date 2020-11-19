#ifndef GUARD_ITERATOR_GROUP_H
#define GUARD_ITERATOR_GROUP_H

#include <iterator>
#include <functional>
#include <type_traits>

namespace mpsym
{

namespace util
{

template<typename T, typename V, bool COPY = false>
class Iterator
{
public:
  using value_type = typename std::remove_const<V>::type;
  using pointer = V *;
  using reference = V &;

private:
  using deref = typename std::conditional<COPY, value_type, reference>::type;

public:
  virtual bool operator==(T const &rhs) const = 0;

  bool operator!=(T const &rhs) const
  { return !((*this) == rhs); }

  T &operator++()
  {
    next();
    return *static_cast<T *>(this);
  }

  T &operator++(int)
  {
    auto ret(*static_cast<T *>(this));
    next();
    return ret;
  }

  deref operator*()
  { return current(); }

  pointer operator->()
  { return &current(); }

private:
  virtual deref current() = 0;
  virtual void next() = 0;
};

template<typename T, typename U>
class IteratorAdaptor
{
  using wrapped = typename T::const_iterator;
  using transform = std::function<U(typename T::const_reference)>;

public:
  class const_iterator : public Iterator<const_iterator, U const, true>
  {
    using base = Iterator<const_iterator, U const>;

  public:
    using value_type = typename base::value_type;
    using pointer = typename base::pointer;
    using reference = typename base::reference;

    const_iterator(wrapped const &it, transform const &f)
    : _it(it),
      _f(f)
    {}

    bool operator==(const_iterator const &rhs) const override
    { return _it == rhs._it; }

  private:
    value_type current() override
    { return _f(*_it); }

    void next() override
    { ++_it; }

    wrapped _it;
    transform _f;
  };

  template<typename FUNC>
  IteratorAdaptor(T &obj, FUNC &&f)
  : _begin(obj.begin()),
    _end(obj.end()),
    _f(f)
  {}

  const_iterator begin()
  { return const_iterator(_begin, _f); }

  const_iterator end()
  { return const_iterator(_end, _f); }

private:
  wrapped _begin, _end;
  transform _f;
};

} // namespace util

} // namespace mpsym

#endif // GUARD_ITERATOR_H
