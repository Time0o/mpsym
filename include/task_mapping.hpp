#ifndef GUARD_TASK_MAPPING_H
#define GUARD_TASK_MAPPING_H

#include <cassert>
#include <initializer_list>
#include <ostream>
#include <type_traits>
#include <utility>
#include <vector>

#include "dump.hpp"
#include "perm.hpp"
#include "util.hpp"

namespace mpsym
{

class TaskMapping : public std::vector<unsigned>
{
public:
  TaskMapping()
  : std::vector<unsigned>()
  {}

  TaskMapping(std::initializer_list<unsigned> tasks)
  : std::vector<unsigned>(tasks)
  {}

  TaskMapping(std::vector<unsigned> tasks)
  : std::vector<unsigned>(tasks)
  {}

  bool less_than(TaskMapping const other) const
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

  template<typename PERM>
  bool less_than(TaskMapping const other,
                 PERM const &perm,
                 unsigned offset = 0u) const
  {
    return foreach_permuted_task(
      perm,
      offset,
      [&](unsigned i, unsigned, unsigned task_permuted, bool &flag){
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

  template<typename PERM>
  void permute(PERM const &perm,
               unsigned offset = 0u,
               bool *modified = nullptr)
  {
    if (modified)
      *modified = false;

    foreach_permuted_task(
      perm,
      offset,
      [&](unsigned i, unsigned task, unsigned task_permuted, bool &){
        if (modified && task_permuted != task)
          *modified = true;

        (*this)[i] = task_permuted;
        return false;
      }
    );
  }

  template<typename PERM>
  TaskMapping permuted(PERM const &perm,
                       unsigned offset = 0u,
                       bool *modified = nullptr) const
  {
    TaskMapping res(*this);

    if (modified)
      *modified = false;

    foreach_permuted_task(
      perm,
      offset,
      [&](unsigned i, unsigned task, unsigned task_permuted, bool &){
        if (modified && task_permuted != task)
          *modified = true;

        res[i] = task_permuted;
        return false;
      }
    );

    return res;
  }

private:
  template<typename PERM, typename FUNC>
  bool foreach_permuted_task_(PERM &&perm,
                              unsigned offset,
                              unsigned degree,
                              FUNC &&func) const
  {
    for (auto i = 0u; i < size(); ++i) {
      unsigned task = (*this)[i];
      if (task <= offset || task > degree + offset)
        continue;

      unsigned task_permuted = perm(task - offset) + offset;

      bool flag;
      if (func(i, task, task_permuted, flag))
        return flag;
    }

    return false;
  }

  template<typename PERM, typename FUNC>
  typename std::enable_if<std::is_same<PERM, internal::Perm>::value, bool>::type
  foreach_permuted_task(PERM const &perm,
                        unsigned offset,
                        FUNC &&func) const
  {
    return foreach_permuted_task_(
      [&](unsigned task){ return perm[task]; },
      offset,
      perm.degree(),
      func);
  }

  template<typename PERM, typename FUNC>
  typename std::enable_if<std::is_same<PERM, internal::PermSet>::value, bool>::type
  foreach_permuted_task(PERM const &perm_word,
                        unsigned offset,
                        FUNC &&func) const
  {
    return foreach_permuted_task_(
      [&](unsigned task){
        for (auto it = perm_word.rbegin(); it != perm_word.rend(); ++it)
          task = (*it)[task];
        return task;
      },
      offset,
      perm_word.degree(),
      func);
  }
};

inline std::ostream &operator<<(std::ostream &os, TaskMapping const &ta)
{
  os << DUMP(static_cast<std::vector<unsigned> const &>(ta));
  return os;
}

} // namespace mpsym

namespace std
{

template<>
struct hash<mpsym::TaskMapping>
{
  std::size_t operator()(mpsym::TaskMapping const &ta) const
  { return mpsym::internal::util::container_hash(ta.begin(), ta.end()); }
};

} // namespace std

#endif // GUARD_TASK_MAPPING_H
