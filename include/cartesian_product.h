#ifndef _GUARD_CARTESIAN_PRODUCT_H
#define _GUARD_CARTESIAN_PRODUCT_H

#include <iterator>

template<typename IT, typename ITBuf>
bool next_in_cartesian_product(IT first, IT last, ITBuf buf_first)
{
  using it_type = typename std::iterator_traits<IT>::value_type::const_iterator;
  using val_type = typename std::iterator_traits<it_type>::value_type;

  auto size = std::distance(first, last);

  static IT current_first, current_last;
  static std::vector<it_type> current_state(size);
  static std::vector<val_type> current_product(size);

  if (first != current_first || last != current_last) {
    current_first = first;
    current_last = last;

    for (auto i = 0u; i < size; ++i) {
      current_state[i] = (first + i)->begin();
      *(buf_first + i) = *current_state[i];
    }

    return true;
  }

  for (auto i = 0u; i < size; ++i) {
    auto it(first + i);

    if (++current_state[i] == it->end())
      current_state[i] = it->begin();

    *(buf_first + i) = *current_state[i];

    if (current_state[i] != it->begin())
      return true;
  }

  return false;
}

#endif // _GUARD_CARTESIAN_PRODUCT_H
