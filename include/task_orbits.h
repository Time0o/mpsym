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

  std::pair<bool, unsigned> insert(TaskMapping const &tm)
  {
    bool new_orbit;
    unsigned equivalence_class;

    auto it = _orbit_index_hash.find(tm.representative);
    if (it == _orbit_index_hash.end()) {
      new_orbit = true;
      equivalence_class = num_orbits();

      _orbit_representatives.push_back(tm.representative);
      _orbit_index_hash[tm.representative] = equivalence_class;

    } else {
      new_orbit = false;
      equivalence_class = it->second;
    }

    return {new_orbit, equivalence_class};
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
