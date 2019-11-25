#ifndef _GUARD_TASK_MAPPING_H
#define _GUARD_TASK_MAPPING_H

#include <ostream>
#include <vector>

#include "util.h"

namespace cgtl
{

using TaskAllocation = std::vector<unsigned>;

struct TaskMappingRequest
{
  TaskMappingRequest(TaskAllocation const &allocation,
                     TaskAllocation::value_type offset = 0u,
                     bool approximate = false)
  : allocation(allocation),
    offset(offset),
    approximate(approximate)
  {}

  TaskAllocation allocation;
  TaskAllocation::value_type offset;
  bool approximate;
};

std::ostream &operator<<(std::ostream &os, TaskMappingRequest const &tmr);

class TaskMapping
{
public:
  TaskMapping(TaskAllocation allocation,
              TaskAllocation representative)
  : allocation(allocation),
    representative(representative)
  {}

  TaskAllocation allocation;
  TaskAllocation representative;
};

std::ostream &operator<<(std::ostream &os, TaskMapping const &tm);

} // namespace cgtl

namespace std
{

template<>
struct hash<cgtl::TaskAllocation>
{
  std::size_t operator()(cgtl::TaskAllocation const &ta) const
  { return util::vector_hash(ta); }
};

} // namespace std

#endif // _GUARD_TASK_MAPPING_H
