#include <algorithm>
#include <fstream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "arch_graph.hpp"
#include "arch_graph_nauty.hpp"
#include "dump.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"

namespace mpsym
{

using namespace internal;

std::string ArchGraph::to_gap() const
{ return graph_nauty().to_gap(num_processors()); }

ArchGraph::ProcessorType ArchGraph::new_processor_type(ProcessorLabel pl)
{
  auto id = _processor_types.size();
  _processor_types.push_back(pl);
  _processor_type_instances.push_back(0u);

  return id;
}

ArchGraph::ChannelType ArchGraph::new_channel_type(ChannelLabel cl)
{
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

void ArchGraph::add_channel(unsigned from, unsigned to, ChannelType cht)
{
  reset_automorphisms();

  _channel_type_instances[cht]++;

  EdgeProperty ep {cht};
  boost::add_edge(from, to, ep, _adj);
}

unsigned ArchGraph::num_processors() const
{ return static_cast<unsigned>(boost::num_vertices(_adj)); }

unsigned ArchGraph::num_channels() const
{ return static_cast<unsigned>(boost::num_edges(_adj)); }

void ArchGraph::dump_processors(std::ostream& os) const
{
  std::vector<std::vector<unsigned>> pes_by_type(_processor_types.size());

  for (auto pe : boost::make_iterator_range(boost::vertices(_adj))) {
    auto pt = _adj[pe].type;

    pes_by_type[pt].push_back(pe);
  }

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
      auto ch = _adj[e].type;

      chs_by_type[ch][pe1].insert(pe2);
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
