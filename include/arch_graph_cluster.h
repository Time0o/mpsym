#ifndef _GUARD_ARCH_GRAPH_CLUSTER_H
#define _GUARD_ARCH_GRAPH_CLUSTER_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "arch_graph_system.h"
#include "bsgs.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"
#include "task_allocation.h"
#include "task_orbits.h"

namespace mpsym
{

class ArchGraphCluster : public ArchGraphSystem
{
public:
  std::string to_gap() const override;

  // TODO: detect equivalent subsystems?
  void add_subsystem(std::shared_ptr<ArchGraphSystem> subsystem)
  { _subsystems.push_back(subsystem); }

  unsigned num_processors() const override;
  unsigned num_channels() const override;
  unsigned num_subsystems() const;

  TaskAllocation mapping(TaskAllocation const &allocation,
                         MappingOptions const *options = nullptr,
                         TaskOrbits *orbits = nullptr) override;

private:
  PermGroup update_automorphisms(BSGS::Options const *bsgs_options) override;

  std::vector<std::shared_ptr<ArchGraphSystem>> _subsystems;
};

} // namespace mpsym

#endif // _GUARD_ARCH_GRAPH_CLUSTER_H
