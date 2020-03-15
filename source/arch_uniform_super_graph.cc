#include <cassert>
#include <memory>
#include <string>

#include "arch_graph_system.h"
#include "arch_uniform_super_graph.h"
#include "bsgs.h"

namespace cgtl
{

std::string
ArchUniformSuperGraph::to_gap() const
{
  return "StandardWreathProduct("
    + _subsystem_super_graph->to_gap() + "," + _subsystem_proto->to_gap() + ")";
}

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
ArchUniformSuperGraph::update_automorphisms(BSGS::Options const *bsgs_options)
{
  return PermGroup::wreath_product(
    _subsystem_proto->automorphisms(bsgs_options),
    _subsystem_super_graph->automorphisms(bsgs_options),
    bsgs_options);
}

} // namespace cgtl
