#ifndef _GUARD_DBG_H
#define _GUARD_DBG_H

#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "dump.h"

class Dbg
{
public:
  enum { TRACE = 1, DBG = 2, INFO = 3, WARN = 4};

  static int loglevel;
  static std::ostream &out;

  Dbg(int level = WARN)
  : _level(level)
  { _buf << _headers[_level]; }

  Dbg &operator<<(char const *val)
  {
    _buf << val;
    return *this;
  }

  template<typename T>
  Dbg &operator<<(T const &val)
  {
    _buf << dump::dump(val);
    return *this;
  }

  ~Dbg() { out << _buf.str() << std::endl; }

private:
  std::string header_indent() {
    return std::string(strlen(_headers[_level]), ' ');
  }

  int _level;
  std::vector<char const *> _headers {
    "", "TRACE: ", "DEBUG: ", "INFO: ", "WARNING: "};
  std::ostringstream _buf;
};

#ifdef NDEBUG

#define Dbg(level) if (0) Dbg()

#else

#define Dbg(level) \
  if (level < Dbg::loglevel) {} \
  else Dbg(level)

#endif

#endif // _GUARD_DBG_H
