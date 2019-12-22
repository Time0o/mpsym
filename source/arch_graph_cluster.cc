#include <cassert>
#include <memory>
#include <vector>

#include "arch_graph_cluster.h"
#include "arch_graph_system.h"
#include "dbg.h"
#include "task_orbits.h"

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

TaskAllocation
ArchGraphCluster::mapping(TaskAllocation const &allocation_,
                          MappingMethod method,
                          TaskOrbits *orbits)
{
  assert(_subsystems.size() > 0u);

  DBG(DEBUG) << "Requested task mapping for: " << allocation_;

  TaskAllocation allocation(allocation_);
  unsigned offset_original = allocation.offset;

  for (auto i = 0u; i < _subsystems.size(); ++i) {
    DBG(DEBUG) << "Subsystem (no. " << i << ")";

    allocation = _subsystems[i]->mapping(allocation, method);

    DBG(DEBUG) << "Yields: " << allocation;

    allocation.offset += _subsystems[i]->num_processors();
  }

  if (orbits)
    orbits->insert(allocation);

  allocation.offset = offset_original;

  return allocation;
}

} // namespace cgtl
