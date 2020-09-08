#ifndef GUARD_PARSE_H
#define GUARD_PARSE_H

#include <string>
#include <vector>

#include "perm.hpp"
#include "perm_set.hpp"
#include "string.hpp"

namespace mpsym
{

namespace util
{

inline internal::Perm parse_perm(unsigned degree, std::string const &str)
{
  std::vector<unsigned> cycle;
  std::vector<std::vector<unsigned>> cycles;

  int n_beg = -1;
  for (int i = 0; i < static_cast<int>(str.size()); ++i) {
    char c = str[i];

    switch (c) {
    case '(':
      cycle.clear();
      break;
    case ',':
    case ')':
      {
        int n = stox<int>(str.substr(n_beg, i - n_beg));

        cycle.push_back(n);
        if (c == ')')
          cycles.push_back(cycle);

        n_beg = -1;
      }
      break;
    default:
      if (n_beg == -1)
        n_beg = i;
    }
  }

  return internal::Perm(degree, cycles);
}

inline internal::PermSet parse_perm_set(unsigned degree,
                                        std::vector<std::string> const &strs)
{
  internal::PermSet ret;
  for (auto const &str : strs)
    ret.insert(parse_perm(degree, str));

  return ret;
}

inline internal::PermSet parse_perm_set(unsigned degree,
                                        std::string const &str)
{
  auto first = str.find('(');
  auto last = str.rfind(')');

  auto str_trimmed(str.substr(first, last - first + 1u));

  auto strs(split(str_trimmed, "),"));

  for (auto &str : strs)
    str += ")";

  return parse_perm_set(degree, strs);
}

} // namespace util

} // namespace mpsym

#endif // GUARD_PARSE_H
