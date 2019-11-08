#ifndef _GUARD_DBG_H
#define _GUARD_DBG_H

#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <vector>

class Dbg
{
public:
  enum { TRACE = 1, DBG = 2, INFO = 3, WARN = 4};

  static int loglevel;
  static std::ostream &out;

  Dbg(int level = WARN) : _level(level) { _buf << _headers[_level]; }

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

  template <typename T>
  Dbg& operator<<(std::unordered_set<T> const &set) {
    if (set.size() == 0u) {
      _buf << "{}";
    } else {
      std::stringstream ss;

      ss << '{';
      for (auto const &x : set)
        ss << x << ", ";
      std::string s(ss.str());
      s.resize(s.size() - 1u);
      s.back() = '}';

      (*this) << s;
    }
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
