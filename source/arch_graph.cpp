#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <nlohmann/json.hpp>

#include "arch_graph.hpp"
#include "dump.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"

using json = nlohmann::json;

namespace mpsym
{

using namespace internal;

std::string ArchGraph::to_gap() const
{ return to_gap_nauty(); }

std::string ArchGraph::to_json() const
{
  // processors dict
  std::map<ProcessorType, std::string> processors_dict;

  // channels dict
  using edge_type = std::pair<ProcessorType, std::string>;
  std::map<ProcessorType, std::vector<edge_type>> channels_dict;

  for (auto pe1 : boost::make_iterator_range(boost::vertices(_adj))) {
    processors_dict[pe1] = _processor_types[_adj[pe1].type];

    for (auto e : boost::make_iterator_range(boost::out_edges(pe1, _adj))) {
      auto pe2 = boost::target(e, _adj);
      channels_dict[pe1].emplace_back(pe2, _channel_types[_adj[e].type]);
    }
  }

  json j_;
  j_["directed"] = _directed;
  j_["processor_types"] = _processor_types;
  j_["channel_types"] = _channel_types;
  j_["processors"] = processors_dict;
  j_["channels"] = channels_dict;

  json j;
  j["graph"] = j_;

  return j.dump();
}

ArchGraph::ProcessorType ArchGraph::new_processor_type(ProcessorLabel pl)
{
  assert(!pl.empty());

  auto id = _processor_types.size();
  _processor_types.push_back(pl);
  _processor_type_instances.push_back(0u);

  return id;
}

ArchGraph::ChannelType ArchGraph::new_channel_type(ChannelLabel cl)
{
  assert(!cl.empty());

  auto id = _channel_types.size();
  _channel_types.push_back(cl);
  _channel_type_instances.push_back(0u);

  return id;
}

unsigned ArchGraph::add_processor(ProcessorType pt)
{
  reset_automorphisms();

  _processor_type_instances[pt]++;

  VertexProperty vp {pt};
  return static_cast<unsigned>(boost::add_vertex(vp, _adj));
}

unsigned ArchGraph::add_processor(ProcessorLabel pl)
{
  ProcessorType pt = assert_processor_type(pl);

  return add_processor(pt);
}

unsigned ArchGraph::add_processors(unsigned n, ProcessorType pt)
{
  assert(n > 0u);

  for (unsigned i = 0u; i < n - 1u; ++i)
    add_processor(pt);

  return add_processor(pt);
}

unsigned ArchGraph::add_processors(unsigned n, ProcessorLabel pl)
{
  assert(n > 0u);

  for (unsigned i = 0u; i < n - 1u; ++i)
    add_processor(pl);

  return add_processor(pl);
}

void ArchGraph::add_channel(unsigned from, unsigned to, ChannelType ct)
{
  if (edge_exists(from, to, ct))
    return;

  reset_automorphisms();

  _channel_type_instances[ct]++;

  EdgeProperty ep {ct};
  boost::add_edge(from, to, ep, _adj);
}

void ArchGraph::add_channel(unsigned pe1, unsigned pe2, ChannelLabel cl)
{
  ChannelType ct = assert_channel_type(cl);

  add_channel(pe1, pe2, ct);
}

void ArchGraph::fully_connect(ChannelType ct)
{
  for (unsigned pe1 = 0u; pe1 < num_processors(); ++pe1) {
    for (unsigned pe2 = (directed() ? 0u : pe1); pe2 < num_processors(); ++pe2)
      add_channel(pe1, pe2, ct);
  }
}

void ArchGraph::fully_connect(ChannelLabel cl)
{
  ChannelType ct = assert_channel_type(cl);

  fully_connect(ct);
}

void ArchGraph::fully_connect(ProcessorType pt, ChannelType ct)
{
  for (unsigned pe1 = 0u; pe1 < num_processors(); ++pe1) {
    if (_adj[pe1].type != pt)
      continue;

    for (unsigned pe2 = pe1 + 1u; pe2 < num_processors(); ++pe2) {
      if (_adj[pe2].type != pt)
        continue;

      add_channel(pe1, pe2, ct);
    }
  }
}

void ArchGraph::fully_connect(ProcessorLabel pl, ChannelLabel cl)
{
  ProcessorType pt = assert_processor_type(pl);
  ChannelType ct = assert_channel_type(cl);

  fully_connect(pt, ct);
}

void ArchGraph::self_connect(ChannelType ct)
{
  for (unsigned pe = 0u; pe < num_processors(); ++pe)
    add_channel(pe, pe, ct);
}

void ArchGraph::self_connect(ChannelLabel cl)
{
  ChannelType ct = assert_channel_type(cl);

  self_connect(ct);
}

void ArchGraph::self_connect(ProcessorType pt, ChannelType ct)
{
  for (unsigned pe = 0u; pe < num_processors(); ++pe) {
    if (_adj[pe].type != pt)
      continue;

    add_channel(pe, pe, ct);
  }
}

void ArchGraph::self_connect(ProcessorLabel pl, ChannelLabel cl)
{
  ProcessorType pt = assert_processor_type(pl);
  ChannelType ct = assert_channel_type(cl);

  self_connect(pt, ct);
}

unsigned ArchGraph::num_processors() const
{ return static_cast<unsigned>(boost::num_vertices(_adj)); }

unsigned ArchGraph::num_channels() const
{ return static_cast<unsigned>(boost::num_edges(_adj)); }

ArchGraph::ChannelType ArchGraph::assert_channel_type(ChannelLabel cl)
{
  ChannelType ct = 0u;
  while (ct < _channel_types.size()) {
    if (_channel_types[ct] == cl)
      break;

    ++ct;
  }

  if (ct == _channel_types.size())
    new_channel_type(cl);

  return ct;
}

ArchGraph::ProcessorType ArchGraph::assert_processor_type(ProcessorLabel pl)
{
  ProcessorType pt = 0u;
  while (pt < _processor_types.size()) {
    if (_processor_types[pt] == pl)
      break;

    ++pt;
  }

  if (pt == _processor_types.size())
    new_processor_type(pl);

  return pt;
}

bool ArchGraph::edge_exists(unsigned from, unsigned to, ChannelType ct) const
{
  return directed() ? edge_exists_directed(from, to, ct)
                    : edge_exists_undirected(from, to, ct);
}

bool ArchGraph::edge_exists_directed(unsigned from, unsigned to, ChannelType ct) const
{
  boost::graph_traits<adjacency_type>::edge_descriptor edge;
  bool edge_exists;

  std::tie(edge, edge_exists) = boost::edge(from, to, _adj);

  return edge_exists && _adj[edge].type == ct;
}

bool ArchGraph::edge_exists_undirected(unsigned from, unsigned to, ChannelType ct) const
{
  return edge_exists_directed(from, to, ct) ||
         edge_exists_directed(to, from, ct);
}

void ArchGraph::dump_processors(std::ostream& os) const
{
  std::vector<std::vector<unsigned>> pes_by_type(_processor_types.size());

  for (auto pe : boost::make_iterator_range(boost::vertices(_adj)))
    pes_by_type[_adj[pe].type].push_back(pe);

  os << "processors: [";

  for (auto pt = 0u; pt < pes_by_type.size(); ++pt) {
    os << "\n  type " << pt;

    auto pt_str = _processor_types[pt];
    if (!pt_str.empty())
      os << " (" << pt_str << ")";

    os << ": " << DUMP(pes_by_type[pt]);
  }

  os << "\n]";
}

void ArchGraph::dump_channels(std::ostream& os) const
{
  std::vector<std::vector<std::set<unsigned>>> chs_by_type(_channel_types.size());

  for (auto &chs : chs_by_type)
    chs.resize(num_processors());

  for (auto pe1 : boost::make_iterator_range(boost::vertices(_adj))) {
    for (auto e : boost::make_iterator_range(boost::out_edges(pe1, _adj))) {
      auto pe2 = boost::target(e, _adj);
      chs_by_type[_adj[e].type][pe1].insert(pe2);
    }
  }

  os << "channels: [";

  for (auto ct = 0u; ct < chs_by_type.size(); ++ct) {
    os << "\n  type " << ct;

    auto ct_str = _channel_types[ct];
    if (!ct_str.empty())
      os << " (" << ct_str << ")";

    os << ": [";

    for (auto pe = 0u; pe < chs_by_type[ct].size(); ++pe) {
      auto adj(chs_by_type[ct][pe]);

      if (adj.empty())
        continue;

      os << "\n    " << pe << ": " << DUMP(adj);
    }

    os << "\n  ]";
  }

  os << "\n]";
}

void ArchGraph::dump_automorphisms(std::ostream& os)
{
  os << "automorphism group: [";

  auto gens(automorphisms().generators());

  for (auto i = 0u; i < gens.size(); ++i) {
    os << "\n  " << gens[i];

    if (i < gens.size() - 1u)
      os << ",";
  }

  os << "\n]";
}

std::ostream &operator<<(std::ostream &os, ArchGraph &ag)
{
  if (ag.num_processors() == 0u) {
    os << "empty architecture graph";
    return os;
  }

  ag.dump_processors(os);
  os << "\n";
  ag.dump_channels(os);
  os << "\n";
  ag.dump_automorphisms(os);
  os << "\n";

  return os;
}

} // namespace mpsym
