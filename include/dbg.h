#ifndef _GUARD_DBG_H
#define _GUARD_DBG_H

#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
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

  ~Dbg()
  { *out << prefix_linebreaks(_buf.str()) << std::endl; }

  Dbg &operator<<(char const *str)
  {
    _buf << str;
    return *this;
  }

  template<typename T>
  Dbg &operator<<(T const &val)
  {
    _buf << DUMP(val);
    return *this;
  }

private:
  std::string header_indent() const
  { return std::string(strlen(_headers[_level]), ' '); }

  std::string prefix_linebreaks(std::string str) const
  {
    std::size_t pos = 0u;

    for (;;) {
      pos = str.find('\n', pos);
      if (pos == std::string::npos)
        break;

      str.replace(pos, 1, "\n" + header_indent());

      ++pos;
    }

    return str;
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
