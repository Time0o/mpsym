#include <cassert>
#include <memory>

#include "arch_graph_system.h"
#include "arch_uniform_super_graph.h"

namespace cgtl
{

unsigned
ArchUniformSuperGraph::num_processors() const
{
  return _subsystem_super_graph->num_processors() *
         _subsystem_proto->num_processors();
}

unsigned
ArchUniformSuperGraph::num_channels() const
{
  return _subsystem_super_graph->num_channels() +
         (_subsystem_super_graph->num_processors() *
          _subsystem_proto->num_channels());
}

PermGroup
ArchUniformSuperGraph::update_automorphisms()
{
  // TODO: swap???
  return PermGroup::wreath_product(
    _subsystem_proto->automorphisms(), _subsystem_super_graph->automorphisms());
}

} // namespace cgtl
