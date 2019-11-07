#include <cassert>
#include <memory>
#include <vector>

#include "arch_graph.h"
#include "dbg.h"

namespace cgtl
{

void
ArchGraphCluster::add_subsystem(std::shared_ptr<ArchGraphSystem> const &ags)
{
  assert(!_automorphisms_valid);

  _subsystems.push_back(ags);
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
    subsystem->complete();

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
ArchGraphCluster::mapping(std::vector<unsigned> const &tasks,
                          unsigned offset,
                          MappingVariant mapping_variant) const
{
  Dbg(Dbg::DBG) << "Requested task mapping for: " << tasks;

  assert(_subsystems.size() > 0u);

  TaskMapping res(tasks, tasks);

  unsigned offs = offset;
  for (auto i = 0u; i < _subsystems.size(); ++i) {
    unsigned next_offs = offs + _subsystems[i]->num_processors();

    Dbg(Dbg::DBG) << "Subsystem (no. " << i << ", "
                  << "pe's " << offs << "-" << next_offs - 1u << ")";

    res = _subsystems[i]->mapping(res.equivalence_class(), offs, mapping_variant);

    Dbg(Dbg::DBG) << "Yields: " << res.equivalence_class();

    offs = next_offs;
  }

  return TaskMapping(tasks, res.equivalence_class());
}

} // namespace cgtl
