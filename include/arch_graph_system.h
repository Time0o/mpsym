#ifndef _GUARD_ARCH_GRAPH_SYSTEM_H
#define _GUARD_ARCH_GRAPH_SYSTEM_H

#include <memory>

#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"
#include "task_allocation.h"
#include "task_orbits.h"

namespace cgtl
{

class ArchGraphSystem
{
public:
  enum class MappingMethod {
    ITERATE,
    LOCAL_SEARCH,
    ORBITS,
    AUTO = ITERATE
  };

  struct MappingOptions {
    MappingMethod method = MappingMethod::AUTO;
    bool match_reprs = true;
  };

  virtual unsigned num_processors() const
  { throw std::logic_error("not implemented"); }

  virtual unsigned num_channels() const
  { throw std::logic_error("not implemented"); }

  virtual PermGroup automorphisms()
  {
    if (!_automorphisms_valid) {
      _automorphisms = update_automorphisms();
      _automorphisms_valid = true;
    }

    return _automorphisms;
  }

  virtual PermSet automorphisms_generators(bool augmented = false);

  virtual PartialPermInverseSemigroup partial_automorphisms()
  { throw std::logic_error("not implemented"); }

  virtual TaskAllocation mapping(TaskAllocation const &allocation,
                                 unsigned offset = 0u,
                                 MappingOptions *options = nullptr,
                                 TaskOrbits *orbits = nullptr);

protected:
  static MappingOptions *get_options(MappingOptions *options)
  { return options ? options : &_default_mapping_options; }

  void set_automorphisms(PermGroup const &automorphisms)
  {
    _automorphisms = automorphisms;
    _automorphisms_valid = true;
  }

  void invalidate_automorphisms()
  { _automorphisms_valid = false; }

private:
  virtual PermGroup update_automorphisms() = 0;

  TaskAllocation min_elem_iterate(TaskAllocation const &tasks,
                                  unsigned offset,
                                  MappingOptions *options,
                                  TaskOrbits *orbits);

  TaskAllocation min_elem_local_search(TaskAllocation const &tasks,
                                       unsigned offset,
                                       MappingOptions *options,
                                       TaskOrbits *orbits);

  TaskAllocation min_elem_orbits(TaskAllocation const &tasks,
                                 unsigned offset,
                                 MappingOptions *options,
                                 TaskOrbits *orbits);

  static bool is_representative(TaskAllocation const &tasks,
                                MappingOptions *options,
                                TaskOrbits *orbits)
  {
    if (!options->match_reprs || !orbits)
      return false;

    return orbits->is_representative(tasks);
  }

  static MappingOptions _default_mapping_options;

  PermGroup _automorphisms;
  bool _automorphisms_valid = false;

  PermSet _augmented_generators;
  bool _augmented_generators_valid = false;
};

// TODO: outstream operator

} // namespace cgtl

#endif // _GUARD_ARCH_GRAPH_SYSTEM_H
