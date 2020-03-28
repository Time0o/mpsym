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
#include "task_allocation.h"
#include "task_orbits.h"

namespace cgtl
{

class ArchUniformSuperGraph : public ArchGraphSystem
{
public:
  ArchUniformSuperGraph(std::shared_ptr<ArchGraphSystem> super_graph,
                        std::shared_ptr<ArchGraphSystem> proto);

  std::string to_gap() const override;

  unsigned num_processors() const override;
  unsigned num_channels() const override;

  TaskAllocation mapping(TaskAllocation const &allocation_,
                         MappingOptions const *options = nullptr,
                         TaskOrbits *orbits = nullptr) override;

private:
  ArchGraphAutomorphisms wreath_product_action_super_graph() const;
  std::vector<ArchGraphAutomorphisms> wreath_product_action_proto() const;

  PermGroup update_automorphisms(BSGS::Options const *bsgs_options) override;

  std::shared_ptr<ArchGraphSystem> _subsystem_super_graph;
  std::shared_ptr<ArchGraphSystem> _subsystem_proto;

  ArchGraphAutomorphisms _sigma_super_graph;
  std::vector<ArchGraphAutomorphisms> _sigmas_proto;
};

} // namespace cgtl

#endif // _GUARD_ARCH_GRAPH_H
