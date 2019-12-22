#ifndef _GUARD_ARCH_GRAPH_CLUSTER_H
#define _GUARD_ARCH_GRAPH_CLUSTER_H

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "arch_graph_system.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"
#include "task_allocation.h"
#include "task_orbits.h"

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

  TaskAllocation mapping(TaskAllocation const &allocation,
                         unsigned offset = 0u,
                         MappingOptions *options = nullptr,
                         TaskOrbits *orbits = nullptr) override;

private:
  PermGroup update_automorphisms() override;

  void add_subsystem(ArchGraphSubsystem &&subsystem);

  std::vector<ArchGraphSubsystem> _subsystems;
};

} // namespace cgtl

#endif // _GUARD_ARCH_GRAPH_CLUSTER_H
