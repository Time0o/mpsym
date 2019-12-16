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
  _automorphisms_valid = false;

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
ArchGraphCluster::update_automorphisms()
{
  assert(!_subsystems.empty());

  _automorphisms = _subsystems[0]->automorphisms();
  for (auto i = 1u; i < _subsystems.size(); ++i) {
    _automorphisms = PermGroup::direct_product(
      _automorphisms, _subsystems[i]->automorphisms());
  }
}

TaskMapping
ArchGraphCluster::mapping(TaskMappingRequest const &tmr_)
{
  assert(_subsystems.size() > 0u);

  DBG(DEBUG) << "Requested task mapping: " << tmr_;

  TaskMappingRequest tmr(tmr_);
  TaskMapping res(tmr.allocation, tmr.allocation);

  for (auto i = 0u; i < _subsystems.size(); ++i) {
    DBG(DEBUG) << "Subsystem (no. " << i << ")";

    auto tm = _subsystems[i]->mapping(tmr);

    DBG(DEBUG) << "Yields: " << tm.representative;

    res.representative = tm.representative;

    tmr.allocation = res.representative;
    tmr.offset += _subsystems[i]->num_processors();
  }

  return res;
}

} // namespace cgtl
