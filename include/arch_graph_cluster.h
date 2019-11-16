#ifndef _GUARD_ARCH_GRAPH_CLUSTER_H
#define _GUARD_ARCH_GRAPH_CLUSTER_H

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "arch_graph_system.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"
#include "task_mapping.h"

namespace cgtl
{

class ArchGraphCluster : public ArchGraphSystem
{
public:
  // TODO: detect equivalent subsystems?
  template<typename ...ARGS>
  void add_subsystem(ARGS &&...args)
  { add_subsystem(ArchGraphSubsystem(std::forward<ARGS>(args)...)); }

  unsigned num_processors() const override;
  unsigned num_channels() const override;

  void complete() override;
  bool completed() const override
  { return _automorphisms_valid; }

  PermGroup automorphisms() const override;
  PartialPermInverseSemigroup partial_automorphisms() const override
  { throw std::logic_error("not implemented"); }

  TaskMapping mapping(TaskMappingRequest const &tmr) const override;

private:
  void add_subsystem(ArchGraphSubsystem &&subsystem);

  PermGroup _automorphisms;
  bool _automorphisms_valid = false;

  std::vector<ArchGraphSubsystem> _subsystems;
};

} // namespace cgtl

#endif // _GUARD_ARCH_GRAPH_CLUSTER_H
