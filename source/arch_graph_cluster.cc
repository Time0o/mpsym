#include <cassert>
#include <memory>
#include <vector>

#include "arch_graph_cluster.h"
#include "arch_graph_system.h"
#include "dbg.h"
#include "task_mapping.h"

namespace cgtl
{

void
ArchGraphCluster::add_subsystem(ArchGraphSubsystem &&subsystem)
{
  assert(!_automorphisms_valid);

  _subsystems.emplace_back(subsystem);
}

unsigned
ArchGraphCluster::num_processors() const
{
  unsigned res = 0u;
  for (auto const &subsystem : _subsystems)
    res += subsystem->num_processors();

  return res;
}

unsigned
ArchGraphCluster::num_channels() const
{
  unsigned res = 0u;
  for (auto const &subsystem : _subsystems)
    res += subsystem->num_channels();

  return res;
}

void
ArchGraphCluster::complete()
{
  assert(!_subsystems.empty());

  for (auto const &subsystem : _subsystems)
    assert(subsystem->completed());

  _automorphisms = _subsystems[0]->automorphisms();
  for (auto i = 1u; i < _subsystems.size(); ++i) {
    _automorphisms = PermGroup::direct_product(
      _automorphisms, _subsystems[i]->automorphisms());
  }

  _automorphisms_valid = true;
}

PermGroup
ArchGraphCluster::automorphisms() const
{
  assert(_automorphisms_valid);

  return _automorphisms;
}

TaskMapping
ArchGraphCluster::mapping(TaskMappingRequest const &tmr_) const
{
  assert(_subsystems.size() > 0u);

  Dbg(Dbg::DBG) << "Requested task mapping: " << tmr_;

  TaskMappingRequest tmr(tmr_);
  TaskMapping res(tmr.allocation, tmr.allocation);

  for (auto i = 0u; i < _subsystems.size(); ++i) {
    Dbg(Dbg::DBG) << "Subsystem (no. " << i << ")";

    auto tm = _subsystems[i]->mapping(tmr);

    Dbg(Dbg::DBG) << "Yields: " << tm.representative;

    res.representative = tm.representative;

    tmr.allocation = res.representative;
    tmr.offset += _subsystems[i]->num_processors();
  }

  return res;
}

} // namespace cgtl
