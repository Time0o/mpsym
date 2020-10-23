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

NautyGraph ArchGraph::graph_nauty() const
{
  int cts = _channel_types.size();
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
    if (level > 0) {
      for (int v = 0; v < n_orig; ++v)
        g.add_edge(v + level * n_orig, v + (level - 1) * n_orig);
    }

    for (auto e : boost::make_iterator_range(boost::edges(_adj))) {
      int source = boost::source(e, _adj);
      int target = boost::target(e, _adj);
      int t = _adj[e].type + 1;

      if (source != target && (t & (1 << level)))
        g.add_edge(source + level * n_orig, target + level * n_orig);
    }
  }

  // set partition
  std::vector<std::vector<int>> ptn(_processor_types.size() * (cts_log2 + 1));

  for (int level = 0; level <= cts_log2; ++level) {
    for (int v = 0; v < n_orig; ++v) {
      int t = _adj[v].type + level * _processor_types.size();
      ptn[t].push_back(v + level * n_orig);
    }
  }

  g.set_partition(ptn);

  return g;
}

std::string ArchGraph::to_gap_nauty() const
{
  auto g(graph_nauty());

  return g.to_gap();
}

PermGroup ArchGraph::automorphisms_nauty(
  AutomorphismOptions const *options) const
{
  auto g(graph_nauty());

  return g.automorphisms(options);
}

} // namespace mpsym
