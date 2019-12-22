#ifndef _GUARD_ARCH_GRAPH_SYSTEM_H
#define _GUARD_ARCH_GRAPH_SYSTEM_H

#include <memory>

#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"
#include "task_allocation.h"
#include "task_orbits.h"

namespace cgtl
{

using ProcessorLabel = char const *;
using ChannelLabel = char const *;
using SubsystemLabel = char const *;

#define DEFAULT_PROCESSOR_LABEL ""
#define DEFAULT_CHANNEL_LABEL ""
#define DEFAULT_SUBSYSTEM_LABEL ""

class ArchGraphSystem
{
public:
  enum class MappingMethod {
    ITERATE,
    LOCAL_SEARCH,
    ORBITS,
    ORBITS_STOP_EARLY,
    AUTO = ITERATE
  };

  struct MappingOptions {
    MappingMethod method = MappingMethod::AUTO;
  };

  ArchGraphSystem()
  : _automorphisms_valid(false),
    _augmented_generators_valid(false)
  {}

  ArchGraphSystem(PermGroup const &automorphisms)
  : _automorphisms(automorphisms),
    _automorphisms_valid(true),
    _augmented_generators_valid(false)
  {}

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
  virtual PermGroup update_automorphisms()
  { throw std::logic_error("not implemented"); }

  TaskAllocation min_elem_iterate(TaskAllocation const &tasks,
                                  unsigned offset);

  TaskAllocation min_elem_local_search(TaskAllocation const &tasks,
                                       unsigned offset);

  TaskAllocation min_elem_orbits(TaskAllocation const &tasks,
                                 unsigned offset,
                                 TaskOrbits *orbits);

  static MappingOptions _default_mapping_options;

  PermGroup _automorphisms;
  bool _automorphisms_valid;

  PermSet _augmented_generators;
  bool _augmented_generators_valid;
};

class ArchGraphSubsystem
{
public:
  ArchGraphSubsystem(std::shared_ptr<ArchGraphSystem> const &ag,
                     SubsystemLabel label = DEFAULT_SUBSYSTEM_LABEL)
  : _arch_graph_system(ag),
    _label(label)
  {}

  template<typename AG>
  ArchGraphSubsystem(AG const &arch_graph_system,
                     SubsystemLabel label = DEFAULT_SUBSYSTEM_LABEL)
  : _arch_graph_system(std::make_shared<AG>(arch_graph_system)),
    _label(label)
  {}

  ArchGraphSystem *operator->() const
  { return _arch_graph_system.get(); }

  SubsystemLabel label() const
  { return _label; }

private:
  std::shared_ptr<ArchGraphSystem> _arch_graph_system;
  SubsystemLabel _label;
};

// TODO: outstream operators for both of these...

} // namespace cgtl

#endif // _GUARD_ARCH_GRAPH_SYSTEM_H
