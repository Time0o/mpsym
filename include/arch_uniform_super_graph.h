#ifndef _GUARD_ARCH_UNIFORM_SUPER_GRAPH_H
#define _GUARD_ARCH_UNIFORM_SUPER_GRAPH_H

#include <memory>
#include <stdexcept>
#include <vector>

#include "arch_graph.h"
#include "arch_graph_system.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"
#include "task_mapping.h"

namespace cgtl
{

class ArchUniformSuperGraph : public ArchGraphSystem
{
public:
  typedef ArchGraph::ProcessorType SubsystemType;
  typedef ArchGraph::ChannelType SubsystemChannelType;

  template<typename ...ARGS>
  ArchUniformSuperGraph(ARGS &&...args)
  : ArchUniformSuperGraph(ArchGraphSubsystem(std::forward<ARGS>(args)...))
  {}

  SubsystemChannelType new_subsystem_channel_type(
    ChannelLabel cl = DEFAULT_CHANNEL_LABEL);

  SubsystemType add_subsystem();
  void add_subsystem_channel(SubsystemType ss1,
                             SubsystemType ss2,
                             SubsystemChannelType ch);

  unsigned num_processors() const override;
  unsigned num_channels() const override;

  void complete() override;
  bool completed() const override
  { return _automorphisms_valid; }

  PermGroup automorphisms() const override;
  PartialPermInverseSemigroup partial_automorphisms() const override
  {  throw std::logic_error("not implemented"); }

  TaskMapping mapping(TaskMappingRequest const &) const override
  { throw std::logic_error("not implemented"); }

private:
  ArchUniformSuperGraph(ArchGraphSubsystem &&subsystem);

  PermGroup _automorphisms;
  bool _automorphisms_valid = false;

  ArchGraph _subsystem_supergraph;
  ArchGraphSubsystem _subsystem_proto;
  ArchGraph::ProcessorType _subsystem_processor_type;
};

} // namespace cgtl

#endif // _GUARD_ARCH_GRAPH_H
