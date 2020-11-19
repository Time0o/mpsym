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
#include "task_mapping_orbit.hpp"

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
  return "FixedPointWreathProduct("
    + _subsystem_proto->to_gap() + ","
    + std::to_string(_subsystem_proto->num_processors()) + ","
    + _subsystem_super_graph->to_gap() + ","
    + std::to_string(_subsystem_super_graph->num_processors()) + ")";
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

std::shared_ptr<ArchGraphAutomorphisms>
ArchUniformSuperGraph::wreath_product_action_super_graph(
  AutomorphismOptions const *options,
  timeout::flag aborted) const
{
  auto automs_super_graph(_subsystem_super_graph->automorphisms());

  unsigned degree_super_graph = _subsystem_super_graph->num_processors();
  unsigned degree_proto = _subsystem_proto->num_processors();
  unsigned degree = degree_super_graph * degree_proto;

  PermSet sigma_super_graph_gens;

  for (auto const &gen_ : automs_super_graph.generators()) {
    std::vector<unsigned> gen(degree);

    for (unsigned i = 0u; i < gen.size(); ++i) {
      unsigned block_from = i / degree_proto;
      unsigned block_offs = i % degree_proto;

      gen[i] = gen_[block_from] * degree_proto + block_offs;
    }

    sigma_super_graph_gens.insert(Perm(gen));
  }

  PermGroup pg(BSGS(degree, sigma_super_graph_gens, options, aborted));

  return std::make_shared<ArchGraphAutomorphisms>(pg);
}

std::vector<std::shared_ptr<ArchGraphAutomorphisms>>
ArchUniformSuperGraph::wreath_product_action_proto(
  AutomorphismOptions const *options,
  timeout::flag aborted) const
{
  auto automs_proto(_subsystem_proto->automorphisms());

  unsigned degree_super_graph = _subsystem_super_graph->num_processors();
  unsigned degree_proto = _subsystem_proto->num_processors();
  unsigned degree = degree_super_graph * degree_proto;

  std::vector<PermSet> sigmas_proto_gens(degree_super_graph);

  for (auto const &gen_ : automs_proto.generators()) {
    for (unsigned b = 0u; b < sigmas_proto_gens.size(); ++b) {
      std::vector<unsigned> gen(degree);
      std::iota(gen.begin(), gen.end(), 0u);

      unsigned block_end = (b + 1u) * degree_proto;
      unsigned block_start = block_end - degree_proto + 1u;

      for (unsigned j = block_start; j <= block_end; ++j)
        gen[j - 1u] = gen_[(j - 1u) % degree_proto] + block_start - 1u;

      sigmas_proto_gens[b].insert(Perm(gen));
    }
  }

  std::vector<std::shared_ptr<ArchGraphAutomorphisms>> sigmas_proto;
  for (auto const &gens : sigmas_proto_gens) {
    PermGroup pg(BSGS(degree, gens, options, aborted));

    sigmas_proto.push_back(std::make_shared<ArchGraphAutomorphisms>(pg));
  }

  return sigmas_proto;
}

PermGroup
ArchUniformSuperGraph::automorphisms_(AutomorphismOptions const *options,
                                      timeout::flag aborted)
{
  return PermGroup::wreath_product(
    _subsystem_proto->automorphisms(options, aborted),
    _subsystem_super_graph->automorphisms(options, aborted),
    options,
    aborted);
}

void
ArchUniformSuperGraph::init_repr_(AutomorphismOptions const *options,
                                  timeout::flag aborted)
{
  auto automs_super_graph(
    _subsystem_super_graph->automorphisms(options, aborted));

  auto automs_proto(
    _subsystem_proto->automorphisms(options, aborted));

  _super_graph_trivial = automs_super_graph.is_trivial();
  _proto_trivial = automs_proto.is_trivial();

  if (_super_graph_trivial || _proto_trivial) {
    _sigma_total = std::make_shared<ArchGraphAutomorphisms>(automorphisms(options, aborted));
  } else {
    _sigma_super_graph = wreath_product_action_super_graph(options, aborted);
    _sigmas_proto = wreath_product_action_proto(options, aborted);
  }

  _sigmas_valid = true;
}

bool
ArchUniformSuperGraph::repr_ready_() const
{
  return _subsystem_super_graph->automorphisms_ready() &&
         _subsystem_proto->automorphisms_ready() &&
         _sigmas_valid;
}

void
ArchUniformSuperGraph::reset_repr_()
{
  _subsystem_super_graph->reset_automorphisms();
  _subsystem_proto->reset_automorphisms();
  _sigmas_valid = false;
}

TaskMapping
ArchUniformSuperGraph::repr_(TaskMapping const &mapping,
                             ReprOptions const *options,
                             TMORs *,
                             timeout::flag aborted)
{
  TaskMapping representative(mapping);

  if (_super_graph_trivial || _proto_trivial)
    return _sigma_total->repr(representative, options, aborted);

  for (auto &sigma : _sigmas_proto)
    representative = sigma->repr(representative, options, aborted);

  return _sigma_super_graph->repr(representative, options, aborted);
}

} // namespace mpsym
