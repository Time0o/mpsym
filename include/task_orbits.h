#ifndef _GUARD_TASK_ORBITS_H
#define _GUARD_TASK_ORBITS_H

#include <algorithm>
#include <cassert>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "task_mapping.h"

namespace mpsym
{

class TaskOrbits
{
  using orbit_reprs_map = std::unordered_map<TaskMapping, unsigned>;

public:
  class const_iterator
  : public std::iterator<std::forward_iterator_tag, TaskMapping>
  {
  public:
    const_iterator(orbit_reprs_map::const_iterator current_it,
                   orbit_reprs_map::const_iterator end_it)
    : _current_it(current_it),
      _end_it(end_it),
      _current_key(current_it == _end_it ? TaskMapping({}) : _current_it->first)
    {}

    const_iterator operator++()
    {
      const_iterator ret = *this;
      next();
      return ret;
    }

    const_iterator operator++(int)
    {
      next();
      return *this;
    }

    TaskMapping const & operator*() const
    { return _current_key; }
    TaskMapping const * operator->() const
    { return &_current_key; }
    bool operator==(const_iterator const &rhs) const
    { return _current_it == rhs._current_it; };
    bool operator!=(const_iterator const &rhs) const
    { return !((*this) == rhs); }

  private:
    void next()
    {
      if (++_current_it != _end_it)
        _current_key = _current_it->first;
    }

    orbit_reprs_map::const_iterator _current_it;
    orbit_reprs_map::const_iterator _end_it;

    TaskMapping _current_key;
  };

  bool operator==(TaskOrbits const &rhs) const
  { return orbit_repr_set() == rhs.orbit_repr_set(); }

  bool operator!=(TaskOrbits const &rhs) const
  { return !(*this == rhs); }

  std::pair<bool, unsigned> insert(TaskMapping const &mapping)
  {
    bool new_orbit;
    unsigned equivalence_class;

    auto it = _orbit_reprs.find(mapping);
    if (it == _orbit_reprs.end()) {
      new_orbit = true;
      equivalence_class = num_orbits();

      _orbit_reprs[mapping] = equivalence_class;

    } else {
      new_orbit = false;
      equivalence_class = it->second;
    }

    return {new_orbit, equivalence_class};
  }

  template<typename IT>
  void insert_all(IT first, IT last)
  {
    for (auto it = first; it != last; ++it)
      insert(*it);
  }

  bool is_repr(TaskMapping const &mapping) const
  {
    auto it(_orbit_reprs.find(mapping));

    return it != _orbit_reprs.end();
  }

  unsigned num_orbits() const
  { return static_cast<unsigned>(_orbit_reprs.size()); }

  const_iterator begin() const
  {
    return const_iterator(_orbit_reprs.begin(),
                          _orbit_reprs.end());
  }

  const_iterator end() const
  {
    return const_iterator(_orbit_reprs.end(),
                          _orbit_reprs.end());
  }

private:
  std::unordered_set<TaskMapping> orbit_repr_set() const
  {
    std::unordered_set<TaskMapping> ret;
    for (auto const &repr : _orbit_reprs)
      ret.insert(repr.first);

    return ret;
  }

  orbit_reprs_map _orbit_reprs;
};

} // namespace mpsym

#endif // _GUARD_TASK_ORBITS_H
