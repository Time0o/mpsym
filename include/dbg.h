#ifndef _GUARD_DBG_H
#define _GUARD_DBG_H

#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "perm.h"
#include "perm_group.h"

class Dbg
{
public:
  enum { TRACE = 1, DBG = 2, INFO = 3, WARN = 4};

  static int loglevel;
  static std::ostream &out;

  Dbg(int level) : _level(level) { _buf << _headers[_level]; }

  template <typename T>
  Dbg& operator<<(T const &val) { _buf << val; return *this; }

  template <typename T>
  Dbg& operator<<(std::vector<T> const &vect) {
    if (vect.size() == 0u) {
      _buf << "[]";
    } else {
      _buf << '[';
      for (auto i = 0u; i < vect.size(); ++i)
        (*this) << vect[i] << ((i == vect.size() - 1u) ? "]" : ", ");
    }
    return *this;
  }

  Dbg& operator<<(cgtl::PermGroup const &pg) {
    _buf << "Permutation Group\n";
    (*this) << header_indent() << "Base: " << pg.bsgs().base << '\n';
    (*this) << header_indent() << "Strong generating set: "
            << pg.bsgs().strong_generators;

    return *this;
  }

  ~Dbg() { _buf << std::endl; out << _buf.str(); }

private:
  std::string header_indent() {
    return std::string(strlen(_headers[_level]), ' ');
  }

  int _level;
  std::vector<char const *> _headers {
    "", "TRACE: ", "DEBUG: ", "INFO: ", "WARNING: "};
  std::ostringstream _buf;
};

#define Dbg(level) \
  if (level < Dbg::loglevel) {} \
  else Dbg(level)

#endif // _GUARD_DBG_H
