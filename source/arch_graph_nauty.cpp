#include <memory>
#include <stdexcept>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "arch_graph.hpp"
#include "arch_graph_nauty.hpp"
#include "arch_graph_system.hpp"
#include "bsgs.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"

namespace mpsym
{

using namespace internal;

PermSet NautyGraph::_gens;
int NautyGraph::_gen_degree;

NautyGraph ArchGraph::graph_nauty() const
{
  int cts = _channel_types.size();
  int cts_log2 = 0; while (cts >>= 1) ++cts_log2;

  int n_orig = num_processors();
  int n = n_orig * (cts_log2 + 1u);

  NautyGraph g(n);

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
      int t = _adj[e].type + 1;

      if (t & (1 << level)){
        int source = boost::source(e, _adj);
        int target = boost::target(e, _adj);

        g.add_edge(source + level * n_orig, target + level * n_orig);
      }
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

PermGroup ArchGraph::automorphisms_nauty(AutomorphismOptions const *options) const
{
  auto g(graph_nauty());
  unsigned g_degree = num_processors();

  if (g_degree == 0u)
    throw std::logic_error("architecture graph has no processing elements");

  return PermGroup(BSGS(num_processors(),
                        g.automorphisms(num_processors()),
                        options));
}

} // namespace mpsym
