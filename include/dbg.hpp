#ifndef GUARD_DBG_H
#define GUARD_DBG_H

#include <cstring>
#include <string>
#include <sstream>
#include <ostream>
#include <vector>

#include "dump.hpp"

#define DBG_NS ::mpsym::internal::dbg

namespace mpsym
{

namespace internal
{

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
  std::stringstream _buf;
};

} // namespace dbg

} // namespace internal

} // namespace mpsym

#ifdef NDEBUG

#define DBG(level) if (0) DBG_NS :: Dbg()

#define DBG_SET_LOGLEVEL(loglevel)

#else

#define DBG(level) \
  if (level < DBG_NS :: Dbg::loglevel) {} \
  else DBG_NS :: Dbg(level)

#define TRACE DBG_NS :: Dbg::TRACE
#define DEBUG DBG_NS :: Dbg::DEBUG
#define INFO DBG_NS :: Dbg::INFO
#define WARN DBG_NS :: Dbg::WARN

#define DBG_SET_LOGLEVEL(level) do { \
  DBG_NS :: Dbg::loglevel = level; } while (0)

#endif

#endif // GUARD_DBG_H
