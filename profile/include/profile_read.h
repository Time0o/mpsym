#ifndef _GUARD_PROFILE_READ_H
#define _GUARD_PROFILE_READ_H

#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>

struct Stream
{
  std::ifstream stream;
  bool valid = false;

  void open(char const *file)
  {
    stream = std::ifstream(file);
    if (!stream)
      throw std::runtime_error("failed to open file");

    valid = true;
  }
};

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

template<typename FUNC>
void foreach_line(std::ifstream &stream, FUNC &&f)
{
  std::string line;
  unsigned lineno = 1;

  while (std::getline(stream, line))
    f(line, lineno++);

  if (stream.bad())
    throw std::runtime_error("failed to read file");
}

#endif // _GUARD_PROFILE_READ_H
