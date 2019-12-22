#ifndef _GUARD_TASK_ALLOCATION_H
#define _GUARD_TASK_ALLOCATION_H

#include <initializer_list>
#include <ostream>
#include <utility>

#include "dump.h"
#include "perm.h"
#include "util.h"

namespace cgtl
{

// TODO: pass min_pe/max_pe during construction
class TaskAllocation : public std::vector<unsigned>
{
public:
  struct SubgraphPerm
  {
    Perm perm;
    unsigned min_pe;
    unsigned max_pe;
  };

  TaskAllocation(std::initializer_list<unsigned> tasks)
  : std::vector<unsigned>(tasks)
  {}

  TaskAllocation(std::vector<unsigned> const &tasks)
  : std::vector<unsigned>(tasks)
  {}

  bool permutes_to_less_than(TaskAllocation const minimum,
                             SubgraphPerm const &perm) const
  {
    return foreach_permuted_task(
      perm,
      [&](unsigned i, unsigned task_permuted, bool &flag){
        unsigned task_min = minimum[i];

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

  void permute(SubgraphPerm const &perm)
  {
    foreach_permuted_task(
      perm,
      [&](unsigned i, unsigned task_permuted, bool &){
        (*this)[i] = task_permuted;
        return false;
      }
    );
  }

  TaskAllocation permuted(SubgraphPerm const &perm) const
  {
    TaskAllocation res(*this);

    foreach_permuted_task(
      perm,
      [&](unsigned i, unsigned task_permuted, bool &){
        res[i] = task_permuted;
        return false;
      }
    );

    return res;
  }

private:
  template<typename FUNC>
  bool foreach_permuted_task(SubgraphPerm const &perm, FUNC &&f) const
  {
    for (auto i = 0u; i < size(); ++i) {
      unsigned task = (*this)[i];
      if (task < perm.min_pe || task > perm.max_pe)
        continue;

      unsigned offs = perm.min_pe - 1u;

      unsigned task_permuted = perm.perm[task - offs] + offs;

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
