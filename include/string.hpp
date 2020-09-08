#ifndef GUARD_STRING_H
#define GUARD_STRING_H

#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
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

template<typename T>
std::string stream(T const &obj)
{
  std::stringstream ss;
  ss << obj;
  return ss.str();
}

template<typename T>
T stox(std::string const &str)
{
  T i;
  bool success;

  try {
    std::size_t idx;

    if (std::is_signed<T>::value)
      i = std::stoll(str, &idx);
    else
      i = std::stoull(str, &idx);

    success = idx == str.size();
  } catch (...) {
    success = false;
  }

  if (!success)
    throw std::invalid_argument("stox failed");

  return i;
}

template<typename T>
T stof(std::string const &str)
{
  T d;
  bool success;

  try {
    std::size_t idx;

    d = std::stod(str, &idx);

    success = idx == str.size();
  } catch (...) {
    success = false;
  }

  if (!success)
    throw std::invalid_argument("stof failed");

  return d;
}

} // namespace util

} // namespace mpsym

#endif // GUARD_STRING_H
