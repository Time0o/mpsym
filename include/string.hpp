#ifndef GUARD_STRING_H
#define GUARD_STRING_H

#include <cstring>
#include <string>
#include <vector>

namespace mpsym
{

namespace util
{

inline std::vector<std::string> split(std::string const &str,
                                      char const *delim = " ")
{
  if (str.find(delim) == std::string::npos)
    return {str};

  std::vector<std::string> res;

  std::size_t pos = 0u, nextpos;

  for (;;) {
    nextpos = str.find(delim, pos);

    if (nextpos == std::string::npos) {
      res.push_back(str.substr(pos, str.size() - pos));
      break;
    } else {
      res.push_back(str.substr(pos, nextpos - pos));
    }

    pos = nextpos + std::strlen(delim);
  }

  return res;
}

inline std::string join(std::vector<std::string> const &strs,
                        char const *delim = ",")
{
  if (strs.empty())
    return "";

  std::string res(strs[0]);
  for (auto i = 1u; i < strs.size(); ++i)
    res += delim + strs[i];

  return res;
}

} // namespace util

} // namespace mpsym

#endif // GUARD_STRING_H
