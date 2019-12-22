#ifndef _GUARD_TASK_MAPPING_H
#define _GUARD_TASK_MAPPING_H

#include <ostream>
#include <vector>

#include "util.h"
#include "task_allocation.h"

namespace cgtl
{

struct TaskMappingRequest
{
  TaskMappingRequest(TaskAllocation const &allocation,
                     TaskAllocation::value_type offset = 0u)
  : allocation(allocation),
    offset(offset)
  {}

  TaskAllocation allocation;
  TaskAllocation::value_type offset;
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

#endif // _GUARD_TASK_MAPPING_H
