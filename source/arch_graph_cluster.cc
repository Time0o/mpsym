#include <cassert>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "arch_graph_cluster.h"
#include "arch_graph_system.h"
#include "bsgs.h"
#include "dbg.h"
#include "task_orbits.h"

namespace mpsym
{

std::string
ArchGraphCluster::to_gap() const
{
  if (_subsystems.empty())
    return "()";

  std::stringstream ss;

  ss << "DirectProduct(" << _subsystems[0];
  for (auto i = 1u; i < _subsystems.size(); ++i)
    ss << "," << _subsystems[i];
  ss << ")";

  return ss.str();
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

unsigned
ArchGraphCluster::num_subsystems() const
{ return static_cast<unsigned>(_subsystems.size()); }


PermGroup
ArchGraphCluster::update_automorphisms(BSGS::Options const *bsgs_options)
{
  assert(!_subsystems.empty());

  std::vector<PermGroup> automorphisms(_subsystems.size());
  for (auto i = 0u; i < _subsystems.size(); ++i)
    automorphisms[i] = _subsystems[i]->automorphisms(bsgs_options);

  return PermGroup::direct_product(automorphisms.begin(),
                                   automorphisms.end(),
                                   bsgs_options);
}

TaskMapping
ArchGraphCluster::repr(TaskMapping const &mapping_,
                       ReprOptions const *options_,
                       TaskOrbits *orbits)
{
  auto options(complete_options(options_));

  assert(_subsystems.size() > 0u);

  DBG(DEBUG) << "Requested task mapping for: " << mapping_;

  TaskMapping mapping(mapping_);

  for (auto i = 0u; i < _subsystems.size(); ++i) {
    DBG(DEBUG) << "Subsystem (no. " << i << ")";

    mapping = _subsystems[i]->repr(mapping, &options);

    DBG(DEBUG) << "Yields: " << mapping;

    options.offset += _subsystems[i]->num_processors();
  }

  if (orbits)
    orbits->insert(mapping);

  return mapping;
}

} // namespace mpsym
