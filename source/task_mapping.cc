#include <algorithm>
#include <ostream>

#include "dump.h"
#include "task_mapping.h"

namespace cgtl
{

std::ostream &operator<<(std::ostream &os, TaskMappingRequest const &tmr)
{
  std::vector<unsigned> ta(tmr.allocation.size());

  std::transform(tmr.allocation.begin(),
                 tmr.allocation.end(),
                 ta.begin(),
                 [&](TaskAllocation::value_type pe){ return pe + tmr.offset; });

  os << DUMP(ta);

  return os;
}

std::ostream &operator<<(std::ostream &os, TaskMapping const &tm)
{
  os << DUMP(tm.allocation) << " => " << DUMP(tm.representative);
  return os;
}

} // namespace cgtl
