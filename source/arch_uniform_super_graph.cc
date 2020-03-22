#include <cassert>
#include <memory>
#include <numeric>
#include <string>

#include "arch_graph_automorphisms.h"
#include "arch_graph_system.h"
#include "arch_uniform_super_graph.h"
#include "block_system.h"
#include "bsgs.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"
#include "task_allocation.h"
#include "task_orbits.h"

namespace cgtl
{

ArchUniformSuperGraph::ArchUniformSuperGraph(
  std::shared_ptr<ArchGraphSystem> super_graph,
  std::shared_ptr<ArchGraphSystem> proto)
: _subsystem_super_graph(super_graph),
  _subsystem_proto(proto),
  _sigma_super_graph(wreath_product_action_super_graph()),
  _sigmas_proto(wreath_product_action_proto())
{}

std::string
ArchUniformSuperGraph::to_gap() const
{
  return "WreathProduct("
    + _subsystem_proto->to_gap() + "," + _subsystem_super_graph->to_gap() + ")";
}

unsigned
ArchUniformSuperGraph::num_processors() const
{
  return _subsystem_super_graph->num_processors() *
         _subsystem_proto->num_processors();
}

unsigned
ArchUniformSuperGraph::num_channels() const
{
  return _subsystem_super_graph->num_channels() +
         (_subsystem_super_graph->num_processors() *
          _subsystem_proto->num_channels());
}

std::vector<ArchGraphAutomorphisms>
ArchUniformSuperGraph::wreath_product_action_proto() const
{
  unsigned degree_super_graph = _subsystem_super_graph->num_processors();
  unsigned degree_proto = _subsystem_proto->num_processors();

  std::vector<PermSet> sigmas_proto_gens(degree_super_graph);

  for (auto const &gen_ : _subsystem_proto->automorphisms_generators()) {
    for (unsigned b = 0u; b < sigmas_proto_gens.size(); ++b) {
      std::vector<unsigned> gen(degree_super_graph * degree_proto);
      std::iota(gen.begin(), gen.end(), 1u);

      unsigned block_end = (b + 1u) * degree_proto;
      unsigned block_start = block_end - degree_proto + 1u;

      for (unsigned j = block_start; j <= block_end; ++j) {
        gen[j - 1u] = gen_[(j - 1u) % degree_proto + 1u] + block_start - 1u;
      }

      sigmas_proto_gens[b].insert(Perm(gen));
    }
  }

  std::vector<ArchGraphAutomorphisms> sigmas_proto;
  for (auto const &gens : sigmas_proto_gens)
    sigmas_proto.emplace_back(PermGroup(gens.degree(), gens));

  return sigmas_proto;
}

ArchGraphAutomorphisms
ArchUniformSuperGraph::wreath_product_action_super_graph() const
{
  unsigned degree_super_graph = _subsystem_super_graph->num_processors();
  unsigned degree_proto = _subsystem_proto->num_processors();

  PermSet sigma_super_graph_gens;

  for (auto const &gen_ : _subsystem_super_graph->automorphisms_generators()) {
    std::vector<unsigned> gen(degree_super_graph * degree_proto);

    for (unsigned i = 1u; i <= gen.size(); ++i) {
      unsigned block_from = (i - 1u) / degree_proto + 1u;
      unsigned block_offs = (i - 1u) % degree_proto;

      unsigned block_to = gen_[block_from];

      gen[i - 1u] = (block_to - 1u) * degree_proto + 1u + block_offs;
    }

    sigma_super_graph_gens.insert(Perm(gen));
  }

  return PermGroup(sigma_super_graph_gens.degree(), sigma_super_graph_gens);
}

PermGroup
ArchUniformSuperGraph::update_automorphisms(BSGS::Options const *bsgs_options)
{
  return PermGroup::wreath_product(
    _subsystem_proto->automorphisms(bsgs_options),
    _subsystem_super_graph->automorphisms(bsgs_options),
    bsgs_options);
}

TaskAllocation
ArchUniformSuperGraph::mapping(TaskAllocation const &allocation_,
                               unsigned offset,
                               MappingOptions *options,
                               TaskOrbits *orbits)
{
  options = get_options(options);

  // mapping
  TaskAllocation allocation(allocation_);

  for (auto &sigma : _sigmas_proto)
    allocation = sigma.mapping(allocation, offset, options);

  allocation = _sigma_super_graph.mapping(allocation, offset, options);

  if (orbits)
    orbits->insert(allocation);

  return allocation;
}

} // namespace cgtl
