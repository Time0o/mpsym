#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iterator>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "arch_graph.h"
#include "arch_graph_system.h"
#include "bsgs.h"
#include "dbg.h"
#include "dump.h"
#include "partial_perm.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"
#include "timer.h"

namespace mpsym
{

std::string ArchGraph::to_gap() const
{
  if (_processor_types.size() > 1u || _channel_types.size() > 1u)
    throw std::logic_error("ArchGraph::to_gap only available for uncolored ArchGraphs");

  typename boost::graph_traits<adjacency_type>::edge_iterator ei, ei_end;

  std::stringstream ss;

  ss << "AutGroupGraph(EdgeOrbitsGraph(Group(()),[";

  for (std::tie(ei, ei_end) = boost::edges(_adj); ei != ei_end; ++ei) {
    auto source = boost::source(*ei, _adj) + 1;
    auto target = boost::target(*ei, _adj) + 1;

    ss << "[" << source << "," << target << "]";

    if (source != target)
      ss << ",[" << target << "," << source << "]";

    if (std::next(ei) != ei_end)
      ss << ",";
  }

  ss << "]," << boost::num_vertices(_adj) << "))";

  return ss.str();
}

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

#if 0
bool ArchGraph::is_partial_automorphism(
  PartialPerm const &pperm) const
{
  for (unsigned x : pperm.dom()) {
    unsigned y = pperm[x];
    if (_adj[x - 1u].type != _adj[y - 1u].type)
      return false;

    for (unsigned z : pperm.dom()) {
      auto edge(boost::edge(x - 1u, z - 1u, _adj));
      if (!std::get<1>(edge))
        continue;

      auto edge_perm(boost::edge(y - 1u, z - 1u, _adj));
      if (!std::get<1>(edge_perm))
        return false;

      auto edge_type = _adj[std::get<0>(edge)].type;
      auto edge_perm_type = _adj[std::get<0>(edge_perm)].type;

      if (edge_type != edge_perm_type)
        return false;
    }
  }

  return true;
}

void ArchGraph::partial_automorphisms_backtrack(
  std::vector<unsigned> dom,
  std::vector<unsigned> im,
  PartialPermInverseSemigroup &res) const
{
  PartialPerm pperm(dom, im);
  DBG(TRACE) << "Considering " << pperm;

  if (!pperm.id()) {
    if (is_partial_automorphism(pperm)) {
      DBG(TRACE) << "=> Is a partial automorphism, adjoining new generator";
      res.adjoin_generators({pperm}, true);
    } else {
      DBG(TRACE) << "=> Is not a partial automorphism, pruning search tree";
      return;
    }
  }

  for (unsigned i = dom.empty() ? 1u : dom.back() + 1u; i <= num_processors(); ++i) {
    dom.push_back(i);

    for (unsigned j = 1u; j <= num_processors(); ++j) {
      if (std::find(im.begin(), im.end(), j) != im.end())
        continue;

      im.push_back(j);

      partial_automorphisms_backtrack(dom, im, res);

      im.pop_back();
    }

    dom.pop_back();
  }
}

PartialPermInverseSemigroup ArchGraph::partial_automorphisms()
{
  DBG(DEBUG) << "=== Determining partial architecture graph automorphisms";

  PartialPermInverseSemigroup res;

  std::vector<unsigned> dom, im;
  dom.reserve(num_processors());
  im.reserve(num_processors());

  partial_automorphisms_backtrack(dom, im, res);

  DBG(DEBUG) << "==> Partial automorphisms are:";
  DBG(DEBUG) << res;

  return res;
}
#endif

ArchGraph ArchGraph::fully_connected(unsigned n,
                                     ProcessorLabel pl,
                                     ChannelLabel cl)
{
  assert(n > 0u);

  ArchGraph ag;

  auto pe(ag.new_processor_type(pl));
  auto ch(ag.new_channel_type(cl));

  std::vector<ArchGraph::ProcessorType> processors;
  for (unsigned i = 0u; i < n; ++i)
    processors.push_back(ag.add_processor(pe));

  for (unsigned i = 0u; i < n; ++i) {
    for (unsigned j = 0u; j < n; ++j) {
      if (j == i)
        continue;

      ag.add_channel(processors[i], processors[j], ch);
    }
  }

  return ag;
}

ArchGraph ArchGraph::regular_mesh(unsigned width,
                                  unsigned height,
                                  ProcessorLabel pl,
                                  ChannelLabel cl)
{
  assert(width > 0u && height > 0u);

  ArchGraph ag;

  auto pe(ag.new_processor_type(pl));
  auto ch(ag.new_channel_type(cl));

  ag.create_mesh(width, height, pe, ch);

  return ag;
}

ArchGraph ArchGraph::hyper_mesh(unsigned width,
                                unsigned height,
                                ProcessorLabel pl,
                                ChannelLabel cl)
{
  assert(width > 0u && height > 0u);

  ArchGraph ag;

  auto pe(ag.new_processor_type(pl));
  auto ch(ag.new_channel_type(cl));

  ag.create_mesh(width, height, pe, ch);

  for (unsigned r = 0u; r < height; ++r) {
    unsigned pe1 = r * width;
    unsigned pe2 = pe1 + width - 1;

    ag.add_channel(pe1, pe2, ch);
  }

  for (unsigned c = 0u; c < width; ++c) {
    unsigned pe1 = c;
    unsigned pe2 = pe1 + (height - 1) * width;

    ag.add_channel(pe1, pe2, ch);
  }

  return ag;
}

void ArchGraph::create_mesh(unsigned width,
                            unsigned height,
                            ProcessorType pe,
                            ChannelType ch)
{
  assert(width > 0u && height > 0u);

  std::vector<std::vector<ArchGraph::ProcessorType>> processors(height);

  for (unsigned i = 0u; i < height; ++i) {
    for (unsigned j = 0u; j < width; ++j) {
      processors[i].push_back(add_processor(pe));
    }
  }

  std::pair<int, int> offsets[] = {{-1, 0}, {0, 1}, {1, 0}, {0, -1}};

  int iheight = static_cast<int>(height);
  int iwidth = static_cast<int>(width);

  for (int i = 0u; i < iheight; ++i) {
    for (int j = 0u; j < iwidth; ++j) {
      for (auto const &offs : offsets) {
        int i_offs = i + offs.first;
        int j_offs = j + offs.second;
        if (i_offs < 0 || i_offs >= iheight || j_offs < 0 || j_offs >= iwidth)
          continue;

        add_channel(processors[i][j], processors[i_offs][j_offs], ch);
      }
    }
  }
}

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
