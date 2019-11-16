#include <cassert>
#include <memory>

#include "arch_graph_system.h"
#include "arch_uniform_super_graph.h"

namespace cgtl
{

ArchUniformSuperGraph::ArchUniformSuperGraph(ArchGraphSubsystem &&subsystem)
: _subsystem_proto(subsystem)
{
  assert(_subsystem_proto->completed());

  _subsystem_processor_type =
     _subsystem_supergraph.new_processor_type(subsystem.label());
}

ArchUniformSuperGraph::SubsystemChannelType
ArchUniformSuperGraph::new_subsystem_channel_type(ChannelLabel cl)
{ return _subsystem_supergraph.new_channel_type(cl); }

ArchUniformSuperGraph::SubsystemType
ArchUniformSuperGraph::add_subsystem()
{
  assert(!_automorphisms_valid);

  return _subsystem_supergraph.add_processor(_subsystem_processor_type);
}

void
ArchUniformSuperGraph::add_subsystem_channel(SubsystemType from,
                                             SubsystemType to,
                                             SubsystemChannelType ch)
{
  assert(!_automorphisms_valid);

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

void
ArchUniformSuperGraph::complete()
{
  assert(_subsystem_proto->completed());

  _subsystem_supergraph.complete();

  _automorphisms = PermGroup::wreath_product(
    _subsystem_proto->automorphisms(), _subsystem_supergraph.automorphisms());

  _automorphisms_valid = true;
}

PermGroup
ArchUniformSuperGraph::automorphisms() const
{
  assert(_automorphisms_valid);

  return _automorphisms;
}

} // namespace cgtl
