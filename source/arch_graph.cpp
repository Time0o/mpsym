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
#include "string.hpp"

using json = nlohmann::json;

namespace mpsym
{

using namespace internal;

std::string ArchGraph::to_gap() const
{ return to_gap_nauty(); }

std::string ArchGraph::to_json() const
{
  // processor_types in use
  decltype(_processor_types) processor_types_in_use;

  for (auto i = 0u; i < num_processor_types(); ++i) {
    if (_processor_type_instances[i] > 0u)
      processor_types_in_use.push_back(_processor_types[i]);
  }

  std::sort(processor_types_in_use.begin(), processor_types_in_use.end());

  // channel types in use
  auto channel_types_in_use(_channel_types);

  std::sort(channel_types_in_use.begin(), channel_types_in_use.end());

  // processors dict
  std::map<ProcessorType, std::string> processors_dict;

  // channels dict
  using edge_type = std::pair<ProcessorType, std::string>;
  std::map<ProcessorType, std::vector<edge_type>> channels_dict;

  for (auto pe : processors()) {
    processors_dict[pe] = processor_type_str(pe);

    auto &channels(channels_dict[pe]);

    for (auto ch : out_channels(pe)) {
      auto channel(std::make_pair(target(ch), channel_type_str(ch)));

      channels.insert(
        std::upper_bound(channels.begin(), channels.end(), channel), channel);
    }
  }

  json j_;
  j_["directed"] = _directed;
  j_["processor_types"] = processor_types_in_use;
  j_["channel_types"] = channel_types_in_use;
  j_["processors"] = processors_dict;
  j_["channels"] = channels_dict;

  json j;
  j["graph"] = j_;

  return j.dump();
}

ArchGraph::ProcessorType ArchGraph::new_processor_type(std::string const &pl)
{
  assert(!pl.empty());

  auto id = num_processor_types();
  _processor_types.push_back(pl);
  _processor_type_instances.push_back(0u);

  return id;
}

ArchGraph::ChannelType ArchGraph::new_channel_type(std::string const &cl)
{
  assert(!cl.empty());

  auto id = num_channel_types();
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

unsigned ArchGraph::add_processor(std::string const &pl)
{
  ProcessorType pt = assert_processor_type(pl);

  return add_processor(pt);
}

void ArchGraph::add_channel(unsigned from, unsigned to, ChannelType ct)
{
  if (channel_exists(from, to, ct))
    return;

  reset_automorphisms();

  _channel_type_instances[ct]++;

  if (from == to)
    add_self_channel(from, ct);
  else
    add_non_self_channel(from, to, ct);
}

void ArchGraph::add_self_channel(unsigned pe, ChannelType ct)
{
  EdgeProperty ep {ct};
  boost::add_edge(pe, pe, ep, _adj);

  auto pt(_adj[pe].type);
  auto pl(_processor_types[pt]);
  auto cl(_channel_types[ct]);

  if (_processor_type_instances[pt] > 0u)
    --_processor_type_instances[pt];

  pt = assert_processor_type(add_self_channel_to_processor_label(pl, cl));

  ++_processor_type_instances[pt];

  _adj[pe].type = pt;
}

void ArchGraph::add_non_self_channel(unsigned from, unsigned to, ChannelType ct)
{
  EdgeProperty ep {ct};
  boost::add_edge(from, to, ep, _adj);
}

std::string ArchGraph::add_self_channel_to_processor_label(
  std::string const &pl, std::string const &cl)
{
  auto tmp(util::split(pl, "%"));

  if (tmp.size() > 1u) {
    auto it(std::lower_bound(tmp.begin() + 1, tmp.end(), cl));
    if (it == tmp.end() || *it != cl)
      tmp.insert(it, cl);
  } else {
    tmp.push_back(cl);
  }

  return util::join(tmp, "%");
}

void ArchGraph::add_channel(unsigned pe1, unsigned pe2, std::string const &cl)
{
  ChannelType ct = assert_channel_type(cl);

  add_channel(pe1, pe2, ct);
}

void ArchGraph::fully_connect(ProcessorType pt, ChannelType ct)
{
  for (unsigned pe1 = 0u; pe1 < num_processors(); ++pe1) {
    if (_adj[pe1].type != pt)
      continue;

    for (unsigned pe2 = (directed() ? 0u : pe1); pe2 < num_processors(); ++pe2) {
      if (_adj[pe2].type != pt)
        continue;

      add_channel(pe1, pe2, ct);
    }
  }
}

void ArchGraph::fully_connect(std::string const &pl, std::string const &cl)
{
  ProcessorType pt = assert_processor_type(pl);
  ChannelType ct = assert_channel_type(cl);

  fully_connect(pt, ct);
}

void ArchGraph::self_connect(ProcessorType pt, ChannelType ct)
{
  for (unsigned pe = 0u; pe < num_processors(); ++pe) {
    if (_adj[pe].type != pt)
      continue;

    add_channel(pe, pe, ct);
  }
}

void ArchGraph::self_connect(std::string const &pl, std::string const &cl)
{
  ProcessorType pt = assert_processor_type(pl);
  ChannelType ct = assert_channel_type(cl);

  self_connect(pt, ct);
}

bool ArchGraph::directed() const
{ return _directed; }

bool ArchGraph::effectively_directed() const
{
  if (!directed())
    return false;

  for (auto ch : channels()) {
    if (!channel_exists(target(ch), source(ch), channel_type(ch)))
      return true;
  }

  return false;
}

unsigned ArchGraph::num_processor_types() const
{ return static_cast<unsigned>(_processor_types.size()); }

unsigned ArchGraph::num_channel_types() const
{ return static_cast<unsigned>(_channel_types.size()); }

unsigned ArchGraph::num_processors() const
{ return static_cast<unsigned>(boost::num_vertices(_adj)); }

unsigned ArchGraph::num_channels() const
{ return static_cast<unsigned>(boost::num_edges(_adj)); }

ArchGraph::ChannelType ArchGraph::assert_channel_type(std::string const &cl)
{
  ChannelType ct = 0u;
  while (ct < num_channel_types()) {
    if (_channel_types[ct] == cl)
      break;

    ++ct;
  }

  if (ct == num_channel_types())
    new_channel_type(cl);

  return ct;
}

ArchGraph::ProcessorType ArchGraph::assert_processor_type(std::string const &pl)
{
  ProcessorType pt = 0u;
  while (pt < num_processor_types()) {
    if (_processor_types[pt] == pl)
      break;

    ++pt;
  }

  if (pt == num_processor_types())
    new_processor_type(pl);

  return pt;
}

bool ArchGraph::channel_exists(unsigned from, unsigned to, ChannelType ct) const
{
  return directed() ? channel_exists_directed(from, to, ct)
                    : channel_exists_undirected(from, to, ct);
}

bool ArchGraph::channel_exists_directed(unsigned from, unsigned to, ChannelType ct) const
{
  for (auto ch : out_channels(from)) {
    if (target(ch) == to && channel_type(ch) == ct)
      return true;
  }

  return false;
}

bool ArchGraph::channel_exists_undirected(unsigned from, unsigned to, ChannelType ct) const
{
  return channel_exists_directed(from, to, ct) ||
         channel_exists_directed(to, from, ct);
}

void ArchGraph::dump_processors(std::ostream& os) const
{
  std::vector<std::vector<unsigned>> pes_by_type(num_processor_types());

  for (auto pe : processors())
    pes_by_type[processor_type(pe)].push_back(pe);

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
  std::vector<std::vector<std::set<unsigned>>> chs_by_type(num_channel_types());

  for (auto &chs : chs_by_type)
    chs.resize(num_processors());

  for (auto pe : processors()) {
    for (auto ch : out_channels(pe))
      chs_by_type[channel_type(ch)][pe].insert(target(ch));
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
