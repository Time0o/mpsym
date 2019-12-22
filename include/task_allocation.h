#ifndef _GUARD_TASK_ALLOCATION_H
#define _GUARD_TASK_ALLOCATION_H

#include <cassert>
#include <initializer_list>
#include <ostream>
#include <utility>
#include <vector>

#include "dump.h"
#include "perm.h"
#include "util.h"

namespace cgtl
{

class TaskAllocation : public std::vector<unsigned>
{
public:
  TaskAllocation(std::initializer_list<unsigned> tasks)
  : std::vector<unsigned>(tasks)
  {}

  TaskAllocation(std::vector<unsigned> const &tasks)
  : std::vector<unsigned>(tasks)
  {}

  bool less_than(TaskAllocation const other) const
  {
    assert(size() == other.size());

    for (auto i = 0u; i < size(); ++i) {
      unsigned task_this = (*this)[i];
      unsigned task_other = other[i];

      if (task_this < task_other)
        return true;
      else if (task_this > task_other)
        return false;
    }

    return false;
  }

  bool less_than(TaskAllocation const other,
                 Perm const &perm,
                 unsigned offset = 0u) const
  {
    return foreach_permuted_task(
      perm,
      offset,
      [&](unsigned i, unsigned task_permuted, bool &flag){
        unsigned task_min = other[i];

        if (task_permuted > task_min) {
          flag = false;
          return true;
        } else if (task_permuted < task_min) {
          flag = true;
          return true;
        }

        return false;
      }
    );
  }

  void permute(Perm const &perm, unsigned offset = 0u)
  {
    foreach_permuted_task(
      perm,
      offset,
      [&](unsigned i, unsigned task_permuted, bool &){
        (*this)[i] = task_permuted;
        return false;
      }
    );
  }

  TaskAllocation permuted(Perm const &perm, unsigned offset = 0u) const
  {
    TaskAllocation res(*this);

    foreach_permuted_task(
      perm,
      offset,
      [&](unsigned i, unsigned task_permuted, bool &){
        res[i] = task_permuted;
        return false;
      }
    );

    return res;
  }

private:
  template<typename FUNC>
  bool foreach_permuted_task(Perm const &perm, unsigned offset, FUNC &&f) const
  {
    for (auto i = 0u; i < size(); ++i) {
      unsigned task = (*this)[i];
      if (task <= offset || task > offset + perm.degree())
        continue;

      unsigned task_permuted = perm[task - offset] + offset;

      bool flag;
      if (f(i, task_permuted, flag))
        return flag;
    }

    return false;
  }
};

inline std::ostream &operator<<(std::ostream &os, TaskAllocation const &ta)
{
  os << DUMP(static_cast<std::vector<unsigned> const &>(ta));
  return os;
}

} // namespace cgtl

namespace std
{

template<>
struct hash<cgtl::TaskAllocation>
{
  std::size_t operator()(cgtl::TaskAllocation const &ta) const
  { return util::container_hash(ta.begin(), ta.end()); }
};

} // namespace std

#endif // _GUARD_TASK_ALLOCATION_H
