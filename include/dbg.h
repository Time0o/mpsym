#ifndef _GUARD_DBG_H
#define _GUARD_DBG_H

#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "dump.h"

namespace dbg
{

class Dbg
{
public:
  enum { TRACE = 1, DEBUG = 2, INFO = 3, WARN = 4};

  static int loglevel;
  static std::ostream *out;

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
    _buf << DUMP(val);
    return *this;
  }

  ~Dbg() { *out << _buf.str() << std::endl; }

private:
  std::string header_indent() {
    return std::string(strlen(_headers[_level]), ' ');
  }

  int _level;
  std::vector<char const *> _headers {
    "", "TRACE: ", "DEBUG: ", "INFO: ", "WARNING: "};
  std::ostringstream _buf;
};

} // namespace dbg

#ifdef NDEBUG

#define DBG(level) if (0) dbg::Dbg()

#define DBG_SET_LOGLEVEL(loglevel)

#define DBG_GET_OUT() dbg::Dbg::out
#define DBG_SET_OUT(os) do { dbg::Dbg::out = os; } while (0)

#else

#define DBG(level) \
  if (level < dbg::Dbg::loglevel) {} \
  else dbg::Dbg(level)

#define TRACE dbg::Dbg::TRACE
#define DEBUG dbg::Dbg::DEBUG
#define INFO dbg::Dbg::INFO
#define WARN dbg::Dbg::WARN

#define DBG_SET_LOGLEVEL(level) do { dbg::Dbg::loglevel = level; } while (0)

#endif

#endif // _GUARD_DBG_H
