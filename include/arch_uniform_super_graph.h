#ifndef _GUARD_ARCH_UNIFORM_SUPER_GRAPH_H
#define _GUARD_ARCH_UNIFORM_SUPER_GRAPH_H

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "arch_graph_automorphisms.h"
#include "arch_graph.h"
#include "arch_graph_system.h"
#include "bsgs.h"
#include "partial_perm_inverse_semigroup.h"
#include "task_mapping.h"
#include "task_orbits.h"

namespace mpsym
{

class ArchUniformSuperGraph : public ArchGraphSystem
{
public:
  ArchUniformSuperGraph(std::shared_ptr<ArchGraphSystem> super_graph,
                        std::shared_ptr<ArchGraphSystem> proto);

  std::string to_gap() const override;

  unsigned num_processors() const override;
  unsigned num_channels() const override;

private:
  PermGroup automorphisms_(AutomorphismOptions const *options) override;

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

  std::shared_ptr<ArchGraphAutomorphisms>
  wreath_product_action_super_graph(AutomorphismOptions const *options) const;

  std::vector<std::shared_ptr<ArchGraphAutomorphisms>>
  wreath_product_action_proto(AutomorphismOptions const *options) const;

  std::shared_ptr<ArchGraphSystem> _subsystem_super_graph;
  std::shared_ptr<ArchGraphSystem> _subsystem_proto;

  std::shared_ptr<ArchGraphAutomorphisms> _sigma_super_graph;
  std::vector<std::shared_ptr<ArchGraphAutomorphisms>> _sigmas_proto;
  bool _sigmas_valid = false;
};

} // namespace mpsym

#endif // _GUARD_ARCH_GRAPH_H
