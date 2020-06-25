#ifndef GUARD_ARCH_GRAPH_CLUSTER_H
#define GUARD_ARCH_GRAPH_CLUSTER_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "arch_graph_system.hpp"
#include "bsgs.hpp"
#include "partial_perm_inverse_semigroup.hpp"
#include "perm_group.hpp"
#include "task_mapping.hpp"
#include "task_orbits.hpp"

namespace mpsym
{

class ArchGraphCluster : public ArchGraphSystem
{
public:
  std::string to_gap() const override;

  // TODO: detect equivalent subsystems?
  void add_subsystem(std::shared_ptr<ArchGraphSystem> subsystem)
  {
    reset_automorphisms();
    _subsystems.push_back(subsystem);
  }

  unsigned num_processors() const override;
  unsigned num_channels() const override;
  unsigned num_subsystems() const;

private:
  internal::BSGS::order_type num_automorphisms_(
    AutomorphismOptions const *options) override
  {
    internal::BSGS::order_type ret = 1;
    for (auto const &subsystem : _subsystems)
      ret *= subsystem->num_automorphisms(options);

    return ret;
  }

  internal::PermGroup automorphisms_(
    AutomorphismOptions const *options) override;

  void init_repr_(AutomorphismOptions const *options) override
  {
    for (auto const &subsystem : _subsystems)
      subsystem->init_repr(options);
  }

  bool repr_ready_() const override
  {
    for (auto const &subsystem : _subsystems) {
      if (!subsystem->repr_ready())
        return false;
    }

    return true;
  }

  void reset_repr_() override
  {
    for (auto const &subsystem : _subsystems)
      subsystem->reset_repr();
  }

  TaskMapping repr_(TaskMapping const &mapping,
                    TaskOrbits *orbits,
                    ReprOptions const *options) override;

  std::vector<std::shared_ptr<ArchGraphSystem>> _subsystems;
};

} // namespace mpsym

#endif // GUARD_ARCH_GRAPH_CLUSTER_H
