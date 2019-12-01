#ifndef _GUARD_PROFILE_READ_H
#define _GUARD_PROFILE_READ_H

#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

inline std::string read_line(std::ifstream &stream)
{
  std::string line;
  if (!std::getline(stream, line))
    throw std::runtime_error("file is empty");

  return line;
}

inline std::string read_file(std::ifstream &stream)
{
  return std::string((std::istreambuf_iterator<char>(stream)),
                      std::istreambuf_iterator<char>());
}

#endif // _GUARD_PROFILE_READ_H
