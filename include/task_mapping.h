#ifndef _GUARD_TASK_MAPPING_H
#define _GUARD_TASK_MAPPING_H

#include <vector>

namespace cgtl
{

class TaskMapping
{
public:
  TaskMapping(std::vector<unsigned> map, std::vector<unsigned> eq)
  : _mapping(map),
    _equivalence_class(eq)
  {}

  std::vector<unsigned> mapping() const
  { return _mapping; }

  std::vector<unsigned> equivalence_class() const
  { return _equivalence_class; }

private:
  std::vector<unsigned> _mapping;
  std::vector<unsigned> _equivalence_class;
};

} // namespace cgtl

#endif // _GUARD_TASK_MAPPING_H
