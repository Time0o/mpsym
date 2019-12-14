#ifndef _GUARD_TASK_ORBITS_H
#define _GUARD_TASK_ORBITS_H

#include <cassert>
#include <unordered_map>
#include <utility>
#include <vector>

#include "task_mapping.h"

namespace cgtl
{

class TaskOrbits
{
public:
  typedef std::vector<TaskAllocation>::const_iterator const_iterator;

  std::pair<bool, unsigned> insert(TaskAllocation const &tr)
  {
    bool new_orbit;
    unsigned equivalence_class;

    auto it = _orbit_index_hash.find(tr);
    if (it == _orbit_index_hash.end()) {
      new_orbit = true;
      equivalence_class = num_orbits();

      _orbit_representatives.push_back(tr);
      _orbit_index_hash[tr] = equivalence_class;

    } else {
      new_orbit = false;
      equivalence_class = it->second;
    }

    return {new_orbit, equivalence_class};
  }

  std::pair<bool, unsigned> insert(TaskMapping const &tm)
  { return insert(tm.representative); }

  template<typename IT>
  void insert_all(IT first, IT last)
  {
    for (auto it = first; it != last; ++it)
      insert(*it);
  }

  unsigned num_orbits() const
  { return static_cast<unsigned>(_orbit_representatives.size()); }

  const_iterator begin() const
  { return _orbit_representatives.begin(); }

  const_iterator end() const
  { return _orbit_representatives.end(); }

private:
  std::vector<TaskAllocation> _orbit_representatives;
  std::unordered_map<TaskAllocation, unsigned> _orbit_index_hash;
};

} // namespace cgtl

#endif // _GUARD_TASK_ORBITS_H
