#ifndef _GUARD_DBG_H
#define _GUARD_DBG_H

#include <iostream>
#include <sstream>
#include <vector>

class Dbg
{
public:
  enum { TRACE = 1, DBG = 2, INFO = 3, WARN = 4};

  static int loglevel;
  static std::ostream &out;

  Dbg(int level) { _buf << _headers[level]; }

  template <typename T>
  Dbg& operator<<(T const &val) { _buf << val; return *this; }

  template <typename T>
  Dbg& operator<<(std::vector<T> const &vect) {
    if (vect.size() == 0u) {
      _buf << "[]";
    } else {
      _buf << '[';
      for (auto i = 0u; i < vect.size(); ++i)
        _buf << vect[i] << ((i == vect.size() - 1u) ? "]" : ", ");
    }
    return *this;
  }

  ~Dbg() { _buf << std::endl; out << _buf.str(); }

private:
  std::vector<char const *> _headers {
    "", "TRACE: ", "DEBUG: ", "INFO: ", "WARNING: "};
  std::ostringstream _buf;
};

#define Dbg(level) \
  if (level < Dbg::loglevel) {} \
  else Dbg(level)

#endif // _GUARD_DBG_H
