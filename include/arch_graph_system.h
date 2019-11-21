#ifndef _GUARD_ARCH_GRAPH_SYSTEM_H
#define _GUARD_ARCH_GRAPH_SYSTEM_H

#include <memory>

#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"
#include "task_mapping.h"

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
  virtual unsigned num_processors() const = 0;
  virtual unsigned num_channels() const = 0;

  virtual PermGroup automorphisms()
  {
    if (!_automorphisms_valid) {
      update_automorphisms();
      _automorphisms_valid = true;
    }

    return _automorphisms;
  }

  virtual PartialPermInverseSemigroup partial_automorphisms()
  { throw std::logic_error("not implemented"); }

  virtual TaskMapping mapping(TaskMappingRequest const &tmr) = 0;

protected:
  PermGroup _automorphisms;
  bool _automorphisms_valid = false;

private:
  virtual void update_automorphisms() = 0;
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
