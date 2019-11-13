#ifndef _GUARD_RESOURCE_UTILITY_H
#define _GUARD_RESOURCE_UTILITY_H

#include <fstream>
#include <stdexcept>
#include <string>

inline std::string resource_path(std::string const &resource)
{
  std::string path("../resources/" + resource);

  std::ifstream is(path.c_str());
  if (!is.good())
    throw std::runtime_error("requested resource does not exist");

  return path;
}

#endif // _GUARD_RESOURCE_UTILITY_H
