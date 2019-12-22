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
  invalidate_automorphisms();

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

PermGroup
ArchGraphCluster::update_automorphisms()
{
  assert(!_subsystems.empty());

  // TODO: optimize
  PermGroup automorphisms(_subsystems[0]->automorphisms());
  for (auto i = 1u; i < _subsystems.size(); ++i)
    automorphisms = PermGroup::direct_product(automorphisms,
                                              _subsystems[i]->automorphisms());

  return automorphisms;
}

TaskAllocation
ArchGraphCluster::mapping(TaskAllocation const &allocation_,
                          unsigned offset,
                          MappingOptions *options,
                          TaskOrbits *orbits)
{
  options = get_options(options);

  assert(_subsystems.size() > 0u);

  DBG(DEBUG) << "Requested task mapping for: " << allocation_;

  TaskAllocation allocation(allocation_);

  for (auto i = 0u; i < _subsystems.size(); ++i) {
    DBG(DEBUG) << "Subsystem (no. " << i << ")";

    allocation = _subsystems[i]->mapping(allocation, offset, options);

    DBG(DEBUG) << "Yields: " << allocation;

    offset += _subsystems[i]->num_processors();
  }

  if (orbits)
    orbits->insert(allocation);

  return allocation;
}

} // namespace cgtl
