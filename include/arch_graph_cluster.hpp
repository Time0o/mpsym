#ifndef GUARD_ARCH_GRAPH_CLUSTER_H
#define GUARD_ARCH_GRAPH_CLUSTER_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "arch_graph_system.hpp"
#include "bsgs.hpp"

namespace mpsym
{

class TaskMapping;
class TaskMapping;

namespace internal { class PermGroup; }

class ArchGraphCluster : public ArchGraphSystem
{
public:
  virtual ~ArchGraphCluster() = default;

  std::string to_gap() const override;
  std::string to_json() const override;

  // TODO: detect equivalent subsystems?
  void add_subsystem(std::shared_ptr<ArchGraphSystem> subsystem)
  {
    reset_automorphisms();
    _subsystems.push_back(subsystem);
  }

  std::vector<std::shared_ptr<ArchGraphSystem>> subsystems() const
  { return _subsystems; }

  unsigned num_processors() const override;
  unsigned num_channels() const override;
  unsigned num_subsystems() const;

private:
  internal::BSGS::order_type num_automorphisms_(
    AutomorphismOptions const *options,
    internal::timeout::aborted_type aborted) override
  {
    internal::BSGS::order_type ret = 1;
    for (auto const &subsystem : _subsystems)
      ret *= subsystem->num_automorphisms(options, aborted);

    return ret;
  }

  internal::PermGroup automorphisms_(AutomorphismOptions const *options,
                                     internal::timeout::aborted_type aborted) override;

  void init_repr_(AutomorphismOptions const *options,
                  internal::timeout::aborted_type aborted) override
  {
    for (auto const &subsystem : _subsystems) {
      if (!subsystem->repr_ready())
        subsystem->init_repr(options, aborted);
    }
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
                    ReprOptions const *options,
                    TaskOrbits *orbits,
                    internal::timeout::aborted_type aborted) override;

  std::vector<std::shared_ptr<ArchGraphSystem>> _subsystems;
};

} // namespace mpsym

#endif // GUARD_ARCH_GRAPH_CLUSTER_H
