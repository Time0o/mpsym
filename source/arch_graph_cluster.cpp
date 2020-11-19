#include <cassert>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "arch_graph_cluster.hpp"
#include "arch_graph_system.hpp"
#include "dump.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"
#include "task_mapping.hpp"
#include "task_mapping_orbit.hpp"

namespace mpsym
{

using namespace internal;

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

std::string
ArchGraphCluster::to_json() const
{
  std::stringstream ss;

  ss << "{\"cluster\": "
     << TRANSFORM_AND_DUMP(_subsystems,
                           [](std::shared_ptr<ArchGraphSystem> const &ags)
                           { return ags->to_json(); })
     << "}";

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
ArchGraphCluster::automorphisms_(AutomorphismOptions const *options,
                                 timeout::flag aborted)
{
  assert(!_subsystems.empty());

  std::vector<PermGroup> automorphisms(_subsystems.size());
  for (auto i = 0u; i < _subsystems.size(); ++i)
    automorphisms[i] = _subsystems[i]->automorphisms(options);

  return PermGroup::direct_product(automorphisms.begin(),
                                   automorphisms.end(),
                                   options);
}

TaskMapping
ArchGraphCluster::repr_(TaskMapping const &mapping_,
                        ReprOptions const *options_,
                        TMORs *,
                        timeout::flag aborted)
{
  auto options(ReprOptions::fill_defaults(options_));

  assert(_subsystems.size() > 0u);

  TaskMapping mapping(mapping_);

  for (auto i = 0u; i < _subsystems.size(); ++i) {
    mapping = _subsystems[i]->repr(mapping, &options, aborted);

    options.offset += _subsystems[i]->num_processors();
  }

  return mapping;
}

} // namespace mpsym
