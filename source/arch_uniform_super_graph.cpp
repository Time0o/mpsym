#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "arch_graph_automorphisms.hpp"
#include "arch_graph_system.hpp"
#include "arch_uniform_super_graph.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"
#include "task_mapping.hpp"
#include "task_orbits.hpp"

namespace mpsym
{

using namespace internal;

ArchUniformSuperGraph::ArchUniformSuperGraph(
  std::shared_ptr<ArchGraphSystem> super_graph,
  std::shared_ptr<ArchGraphSystem> proto)
: _subsystem_super_graph(super_graph),
  _subsystem_proto(proto)
{}

std::string
ArchUniformSuperGraph::to_gap() const
{
  return "WreathProduct("
    + _subsystem_proto->to_gap() + "," + _subsystem_super_graph->to_gap() + ")";
}

std::string
ArchUniformSuperGraph::to_json() const
{
  std::stringstream ss;

  ss << "{\"super_graph\": ["
     << _subsystem_proto->to_json()
     << ", "
     << _subsystem_super_graph->to_json()
     << "]}";

  return ss.str();
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
  auto inter_channels =
    (_subsystem_proto->num_processors() * _subsystem_proto->num_processors()) *
     _subsystem_super_graph->num_channels();

  auto intra_channels = _subsystem_super_graph->num_processors() *
                        _subsystem_proto->num_channels();

  return inter_channels + intra_channels;
}

std::vector<std::shared_ptr<ArchGraphAutomorphisms>>
ArchUniformSuperGraph::wreath_product_action_proto(
  AutomorphismOptions const *options) const
{
  unsigned degree_super_graph = _subsystem_super_graph->num_processors();
  unsigned degree_proto = _subsystem_proto->num_processors();

  std::vector<PermSet> sigmas_proto_gens(degree_super_graph);

  for (auto const &gen_ :
       _subsystem_proto->automorphisms(options).generators()) {

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

  std::vector<std::shared_ptr<ArchGraphAutomorphisms>> sigmas_proto;

  for (auto const &gens : sigmas_proto_gens) {
    PermGroup pg(gens.degree(), gens);

    sigmas_proto.push_back(std::make_shared<ArchGraphAutomorphisms>(pg));
  }

  return sigmas_proto;
}

std::shared_ptr<ArchGraphAutomorphisms>
ArchUniformSuperGraph::wreath_product_action_super_graph(
  AutomorphismOptions const *options) const
{
  unsigned degree_super_graph = _subsystem_super_graph->num_processors();
  unsigned degree_proto = _subsystem_proto->num_processors();

  PermSet sigma_super_graph_gens;

  for (auto const &gen_ :
       _subsystem_super_graph->automorphisms(options).generators()) {

    std::vector<unsigned> gen(degree_super_graph * degree_proto);

    for (unsigned i = 1u; i <= gen.size(); ++i) {
      unsigned block_from = (i - 1u) / degree_proto + 1u;
      unsigned block_offs = (i - 1u) % degree_proto;

      unsigned block_to = gen_[block_from];

      gen[i - 1u] = (block_to - 1u) * degree_proto + 1u + block_offs;
    }

    sigma_super_graph_gens.insert(Perm(gen));
  }

  PermGroup pg(sigma_super_graph_gens.degree(), sigma_super_graph_gens);

  return std::make_shared<ArchGraphAutomorphisms>(pg);
}

PermGroup
ArchUniformSuperGraph::automorphisms_(AutomorphismOptions const *options)
{
  return PermGroup::wreath_product(
    _subsystem_proto->automorphisms(options),
    _subsystem_super_graph->automorphisms(options),
    options);
}

TaskMapping
ArchUniformSuperGraph::repr_(TaskMapping const &mapping_,
                             ReprOptions const *options,
                             TaskOrbits *)
{
  TaskMapping mapping(mapping_);

  for (auto &sigma : _sigmas_proto)
    mapping = sigma->repr(mapping, options);

  mapping = _sigma_super_graph->repr(mapping, options);

  return mapping;
}

} // namespace mpsym
