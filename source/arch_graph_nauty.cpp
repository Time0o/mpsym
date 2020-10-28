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
  std::vector<std::vector<int>> ptn(num_processor_types() * (cts_log2 + 1));

  for (int level = 0; level <= cts_log2; ++level) {
    for (int v = 0; v < n_orig; ++v)
      ptn[processor_type(v) + level * num_processor_types()].push_back(v + level * n_orig);
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
