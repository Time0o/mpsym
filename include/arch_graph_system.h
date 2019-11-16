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

  virtual void complete() = 0;
  virtual bool completed() const = 0;

  virtual PermGroup automorphisms() const = 0;
  virtual PartialPermInverseSemigroup partial_automorphisms() const = 0;

  virtual TaskMapping mapping(TaskMappingRequest const &tmr) const = 0;
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
