#include <cassert>
#include <memory>
#include <string>

#include "arch_graph.h"

namespace cgtl
{

ArchUniformSuperGraph::ArchUniformSuperGraph(
  std::shared_ptr<ArchGraphSystem> const &subsystem_proto,
  std::string const &subsystem_label)
: _subsystem_proto(subsystem_proto)
{
  _subsystem_proto->complete();

  _subsystem_processor_type =
     _subsystem_supergraph.new_processor_type(subsystem_label);
}

void
ArchUniformSuperGraph::set_subsystem(
  std::shared_ptr<ArchGraph> const &subsystem_proto,
  std::string const &subsystem_label)
{
  _subsystem_proto = subsystem_proto;

  _subsystem_proto->complete();

  _subsystem_processor_type =
     _subsystem_supergraph.new_processor_type(subsystem_label);
}

ArchUniformSuperGraph::SubsystemChannelType
ArchUniformSuperGraph::new_subsystem_channel_type(std::string const &label)
{
  return _subsystem_supergraph.new_channel_type(label);
}

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
