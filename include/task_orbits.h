#ifndef _GUARD_TASK_ORBITS_H
#define _GUARD_TASK_ORBITS_H

#include <algorithm>
#include <cassert>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "task_allocation.h"

namespace cgtl
{

class TaskOrbits
{
  using orbit_representatives_map = std::unordered_map<TaskAllocation, unsigned>;

public:
  class const_iterator
  : public std::iterator<std::forward_iterator_tag, TaskAllocation>
  {
  public:
    const_iterator(orbit_representatives_map::const_iterator current_it,
                   orbit_representatives_map::const_iterator end_it)
    : _current_it(current_it),
      _end_it(end_it),
      _current_key(current_it == _end_it ? TaskAllocation({}) : _current_it->first)
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

    TaskAllocation const & operator*() const
    { return _current_key; }
    TaskAllocation const * operator->() const
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

    orbit_representatives_map::const_iterator _current_it;
    orbit_representatives_map::const_iterator _end_it;

    TaskAllocation _current_key;
  };

  bool operator==(TaskOrbits const &rhs) const
  { return orbit_representative_set() == rhs.orbit_representative_set(); }

  bool operator!=(TaskOrbits const &rhs) const
  { return !(*this == rhs); }

  std::pair<bool, unsigned> insert(TaskAllocation const &allocation)
  {
    bool new_orbit;
    unsigned equivalence_class;

    auto it = _orbit_representatives.find(allocation);
    if (it == _orbit_representatives.end()) {
      new_orbit = true;
      equivalence_class = num_orbits();

      _orbit_representatives[allocation] = equivalence_class;

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

  bool is_representative(TaskAllocation const &allocation) const
  { return std::find(begin(), end(), allocation) != end(); }

  unsigned num_orbits() const
  { return static_cast<unsigned>(_orbit_representatives.size()); }

  const_iterator begin() const
  {
    return const_iterator(_orbit_representatives.begin(),
                          _orbit_representatives.end());
  }

  const_iterator end() const
  {
    return const_iterator(_orbit_representatives.end(),
                          _orbit_representatives.end());
  }

private:
  std::unordered_set<TaskAllocation> orbit_representative_set() const
  {
    std::unordered_set<TaskAllocation> ret;
    for (auto const &repr : _orbit_representatives)
      ret.insert(repr.first);

    return ret;
  }

  orbit_representatives_map _orbit_representatives;
};

} // namespace cgtl

#endif // _GUARD_TASK_ORBITS_H
