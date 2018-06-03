#ifndef _GUARD_DBG_H
#define _GUARD_DBG_H

#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>

#include "perm.h"

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
    (*this) << header_indent() << "Base: "
            << pg.base() << '\n';
    (*this) << header_indent() << "Strong generating set: "
            << pg.strong_generating_set() << '\n';

    std::vector<std::vector<unsigned>> orbits;
    std::vector<std::vector<cgtl::Perm>> transversals;
    for (auto const &st : pg.schreier_trees()) {
      orbits.push_back(st.orbit());
      transversals.push_back(st.transversals(orbits.back()));
    }

    (*this) << header_indent() << "Orbits: " << orbits << '\n';
    (*this) << header_indent() << "Transversals: " << transversals;

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
