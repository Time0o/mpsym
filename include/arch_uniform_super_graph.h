#ifndef _GUARD_ARCH_UNIFORM_SUPER_GRAPH_H
#define _GUARD_ARCH_UNIFORM_SUPER_GRAPH_H

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "arch_graph.h"
#include "arch_graph_system.h"
#include "bsgs.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"

namespace cgtl
{

class ArchUniformSuperGraph : public ArchGraphSystem
{
public:
  ArchUniformSuperGraph(std::shared_ptr<ArchGraphSystem> super_graph,
                        std::shared_ptr<ArchGraphSystem> proto)
  : _subsystem_super_graph(super_graph),
    _subsystem_proto(proto)
  {}

  std::string to_gap() const override;

  unsigned num_processors() const override;
  unsigned num_channels() const override;

private:
  PermGroup update_automorphisms(BSGS::Options const *bsgs_options) override;

  std::shared_ptr<ArchGraphSystem> _subsystem_super_graph;
  std::shared_ptr<ArchGraphSystem> _subsystem_proto;
};

} // namespace cgtl

#endif // _GUARD_ARCH_GRAPH_H
