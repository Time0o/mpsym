#ifndef GUARD_ITERATOR_GROUP_H
#define GUARD_ITERATOR_GROUP_H

#include <iterator>
#include <type_traits>

namespace mpsym
{

namespace internal
{

template<typename T, typename V>
class ForwardIterator
{
public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = typename std::remove_const<V>::type;
  using pointer = V *;
  using reference = V &;

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

  reference operator*()
  { return current(); }

  pointer operator->()
  { return &current(); }

private:
  virtual reference current() = 0;
  virtual void next() = 0;
};

} // namespace internal

} // namespace mpsym

#endif // GUARD_ITERATOR_H
