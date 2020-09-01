#ifndef GUARD_ARCH_UNIFORM_SUPER_GRAPH_H
#define GUARD_ARCH_UNIFORM_SUPER_GRAPH_H

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "arch_graph_system.hpp"
#include "bsgs.hpp"

namespace mpsym
{

class TaskMapping;
class TaskOrbits;

namespace internal { class ArchGraphAutomorphisms; }

class ArchUniformSuperGraph : public ArchGraphSystem
{
public:
  virtual ~ArchUniformSuperGraph() = default;

  ArchUniformSuperGraph(std::shared_ptr<ArchGraphSystem> super_graph,
                        std::shared_ptr<ArchGraphSystem> proto);

  std::string to_gap() const override;

  std::string to_json() override;

  unsigned num_processors() const override;
  unsigned num_channels() const override;

private:
  internal::BSGS::order_type num_automorphisms_(
    AutomorphismOptions const *options) override
  {
    using boost::multiprecision::pow;

    auto order_super_graph(_subsystem_super_graph->num_automorphisms(options));
    auto order_proto(_subsystem_proto->num_automorphisms(options));

    unsigned lmp_super_graph =
      _subsystem_super_graph->automorphisms().generators().largest_moved_point();

    return pow(order_proto, lmp_super_graph) * order_super_graph;
  }

  internal::PermGroup automorphisms_(
    AutomorphismOptions const *options) override;

  void init_repr_(AutomorphismOptions const *options) override
  {
    _sigma_super_graph = wreath_product_action_super_graph(options);
    _sigmas_proto = wreath_product_action_proto(options);
    _sigmas_valid = true;
  }

  bool repr_ready_() const override
  {
    return _subsystem_super_graph->automorphisms_ready() &&
           _subsystem_proto->automorphisms_ready() &&
           _sigmas_valid;
  }

  void reset_repr_() override
  {
    _subsystem_super_graph->reset_automorphisms();
    _subsystem_proto->reset_automorphisms();
    _sigmas_valid = false;
  }

  TaskMapping repr_(TaskMapping const &mapping_,
                    TaskOrbits *orbits,
                    ReprOptions const *options) override;

  std::shared_ptr<internal::ArchGraphAutomorphisms>
  wreath_product_action_super_graph(AutomorphismOptions const *options) const;

  std::vector<std::shared_ptr<internal::ArchGraphAutomorphisms>>
  wreath_product_action_proto(AutomorphismOptions const *options) const;

  std::shared_ptr<ArchGraphSystem> _subsystem_super_graph;
  std::shared_ptr<ArchGraphSystem> _subsystem_proto;

  std::shared_ptr<internal::ArchGraphAutomorphisms> _sigma_super_graph;
  std::vector<std::shared_ptr<internal::ArchGraphAutomorphisms>> _sigmas_proto;
  bool _sigmas_valid = false;
};

} // namespace mpsym

#endif // GUARD_ARCH_GRAPH_H
