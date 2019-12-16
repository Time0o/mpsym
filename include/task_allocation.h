#include <initializer_list>
#include <ostream>
#include <utility>

#include "dump.h"
#include "perm.h"
#include "util.h"

namespace cgtl
{

// TODO: pass min_pe/max_pe during construction
struct TaskAllocation : public std::vector<unsigned>
{
  TaskAllocation(std::initializer_list<unsigned> tasks)
  : std::vector<unsigned>(tasks)
  {}

  TaskAllocation(std::vector<unsigned> const &tasks)
  : std::vector<unsigned>(tasks)
  {}

  bool minimizes(Perm const &perm, unsigned min_pe, unsigned max_pe) const
  {
    for (unsigned task : *this) {
      if (task < min_pe || task > max_pe)
        continue;

      unsigned offs = min_pe - 1u;

      unsigned permuted = perm[task - offs] + offs;

      if (permuted > task)
        return false;
      else if (permuted < task)
        return true;
    }

    return false;
  }

  void permute(Perm const &perm, unsigned min_pe, unsigned max_pe)
  {
    for (unsigned &task : *this) {
      if (task < min_pe || task > max_pe)
        continue;

      unsigned offs = min_pe - 1u;

      task = perm[task - offs] + offs;
    }
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
