#include <cassert>
#include <memory>

#include "arch_graph_system.h"
#include "arch_uniform_super_graph.h"

namespace cgtl
{

ArchUniformSuperGraph::ArchUniformSuperGraph(ArchGraphSubsystem &&subsystem)
: _subsystem_proto(subsystem)
{
  _subsystem_processor_type =
     _subsystem_supergraph.new_processor_type(subsystem.label());
}

ArchUniformSuperGraph::SubsystemChannelType
ArchUniformSuperGraph::new_subsystem_channel_type(ChannelLabel cl)
{ return _subsystem_supergraph.new_channel_type(cl); }

ArchUniformSuperGraph::SubsystemType
ArchUniformSuperGraph::add_subsystem()
{
  invalidate_automorphisms();

  return _subsystem_supergraph.add_processor(_subsystem_processor_type);
}

void
ArchUniformSuperGraph::add_subsystem_channel(SubsystemType from,
                                             SubsystemType to,
                                             SubsystemChannelType ch)
{
  invalidate_automorphisms();

  _subsystem_supergraph.add_channel(from, to, ch);
}

unsigned
ArchUniformSuperGraph::num_processors() const
{
  return _subsystem_supergraph.num_processors() *
         _subsystem_proto->num_processors();
}

unsigned
ArchUniformSuperGraph::num_channels() const
{
  unsigned res = _subsystem_supergraph.num_channels();

  res += _subsystem_supergraph.num_processors() *
         _subsystem_proto->num_channels();

  return res;
}

PermGroup
ArchUniformSuperGraph::update_automorphisms()
{
  return PermGroup::wreath_product(
    _subsystem_proto->automorphisms(), _subsystem_supergraph.automorphisms());
}

} // namespace cgtl
