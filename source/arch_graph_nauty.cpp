#include <iostream>
#include <string>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>

extern "C" {
  #include "nauty.h"
}

#include "arch_graph.hpp"
#include "nauty_graph.hpp"
#include "perm_group.hpp"

namespace mpsym
{

using namespace internal;

NautyGraph ArchGraph::graph_nauty(bool resolve_loops) const
{
  int cts = num_channel_types();
  int cts_log2 = 0; while (cts >>= 1) ++cts_log2;

  int n_orig = num_processors();
  int n = n_orig * (cts_log2 + 1u);

  NautyGraph g(n, n_orig, directed(), effectively_directed());

  /* node numbering:
   *  ...     ...           ...
   *   |       |             |
   * (n+1)---(n+2)-- ... --(n+n)
   *   |       |             |
   *  (1)-----(2)--- ... ---(n)
   */

  // add edges
  for (int level = 0; level <= cts_log2; ++level) {
    for (int level_ = 0; level_ <= cts_log2; ++level_) {
      if (level_ == level)
        continue;

      for (int v = 0; v < n_orig; ++v)
        g.add_edge(v + level * n_orig, v + level_ * n_orig);
    }

    for (auto ch : channels()) {
      if (source(ch) != target(ch) && (channel_type(ch) + 1 & (1 << level)))
        g.add_edge(source(ch) + level * n_orig, target(ch) + level * n_orig);
    }
  }

  // set partition
  std::vector<std::vector<int>> ptn(
    num_processor_types(resolve_loops) * (cts_log2 + 1));

  for (int level = 0; level <= cts_log2; ++level) {
    for (int v = 0; v < n_orig; ++v) {
      auto p = processor_type(v, resolve_loops)
               + level * num_processor_types(resolve_loops);

      ptn[p].push_back(v + level * n_orig);
    }
  }

  g.set_partition(ptn);

  return g;
}

void ArchGraph::resolve_loops_nauty()
{
  std::map<std::string, std::vector<unsigned>> processor_types_no_loops;

  for (auto pe : processors()) {
    // append loop channel names to processor name
    std::vector<std::string> dummy_pt_components{processor_type_str(pe)};

    for (auto ch : out_channels(pe)) {
      if (source(ch) != target(ch))
        continue;

      dummy_pt_components.push_back(channel_type_str(ch));
    }

    if (dummy_pt_components.size() > 1u)
      std::sort(dummy_pt_components.begin() + 1, dummy_pt_components.end());

    auto dummy_pt(util::join(dummy_pt_components, "%"));

    processor_types_no_loops[dummy_pt].push_back(pe);
  }

  ProcessorType pt = 0u;
  for (auto const &tmp : processor_types_no_loops) {
    for (unsigned pe : tmp.second)
      _processor_types_loops_resolved[pe] = pt;

    ++pt;
  }
}

std::string ArchGraph::to_gap_nauty() const
{
  auto g(graph_nauty(false));

  return g.to_gap();
}

PermGroup ArchGraph::automorphisms_nauty(
  AutomorphismOptions const *options) const
{
  auto g(graph_nauty(true));

  return g.automorphisms(options);
}

} // namespace mpsym
