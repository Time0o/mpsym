#include <vector>

#include "arch_graph.h"
#include "bsgs.h"
#include "nauty.h"
#include "perm_set.h"

namespace
{

int nauty_generator_degree;
mpsym::PermSet nauty_generators;

void nauty_save_generator(int, int *perm, int *, int, int, int)
{
  std::vector<unsigned> tmp(nauty_generator_degree);
  for (int i = 0; i < nauty_generator_degree; ++i)
    tmp[i] = perm[i] + 1;

  nauty_generators.emplace(tmp);
}

void nauty_free()
{
  naugraph_freedyn();
  nautil_freedyn();
  nauty_freedyn();
}

} // anonymous namespace

namespace mpsym
{

PermGroup ArchGraph::automorphisms_nauty(AutomorphismOptions const *options)
{
  // allocate nauty structures
  DYNALLSTAT(graph, g, g_sz);
  DYNALLSTAT(int, lab, lab_sz);
  DYNALLSTAT(int, ptn, ptn_sz);
  DYNALLSTAT(int, orbits, orbits_sz);

  static DEFAULTOPTIONS_GRAPH(nauty_options);
  nauty_options.defaultptn = FALSE;
  nauty_options.userautomproc = nauty_save_generator;

  statsblk stats;

  int cts = static_cast<int>(_channel_types.size());
  int cts_log2 = 0; while (cts >>= 1) ++cts_log2;

  int n_orig = static_cast<int>(boost::num_vertices(_adj));
  int n = static_cast<int>(n_orig * (cts_log2 + 1u));
  int m = SETWORDSNEEDED(n);

#ifndef NDEBUG
  nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
#endif

  DYNALLOC2(graph, g, g_sz, m, n, "malloc graph");
  DYNALLOC1(int, lab, lab_sz, n, "malloc lab");
  DYNALLOC1(int, ptn, ptn_sz, n, "malloc ptn");
  DYNALLOC1(int, orbits, orbits_sz, n, "malloc orbits");

  // construct nauty graph
  EMPTYGRAPH(g, m, n);

  /* node numbering:
   *  ...     ...           ...
   *   |       |             |
   * (n+1)---(n+2)-- ... --(n+n)
   *   |       |             |
   *  (1)-----(2)--- ... ---(n)
   */
  for (int level = 0; level <= cts_log2; ++level) {
    if (level > 0) {
      for (int v = 0; v < n_orig; ++v)
        ADDONEEDGE(g, v + level * n_orig, v + (level - 1) * n_orig, m);
    }

    for (auto e : boost::make_iterator_range(boost::edges(_adj))) {
      int t = static_cast<int>(_adj[e].type) + 1;

      if (t & (1 << level)){
        int source = static_cast<int>(boost::source(e, _adj));
        int target = static_cast<int>(boost::target(e, _adj));

        ADDONEEDGE(g, source + level * n_orig, target + level * n_orig, m);
      }
    }
  }

  std::vector<std::vector<int>> ptn_(_processor_types.size() * (cts_log2 + 1));

  for (int level = 0; level <= cts_log2; ++level) {
    for (int v = 0; v < n_orig; ++v) {
      int t = static_cast<int>(_adj[v].type + level * _processor_types.size());
      ptn_[t].push_back(v + level * n_orig);
    }
  }

  int i = 0;
  for (auto const &p : ptn_) {
    for (auto j = 0u; j < p.size(); ++j) {
      lab[i] = p[j];
      ptn[i] = (j == p.size() - 1u) ? 0 : 1;
      ++i;
    }
  }

  // call nauty
  nauty_generators.clear();
  nauty_generator_degree = n_orig;

  densenauty(g, lab, ptn, orbits, &nauty_options, &stats, m, n, nullptr);

  DYNFREE(g, g_sz);
  DYNFREE(lab, lab_sz);
  DYNFREE(ptn, ptn_sz);
  DYNFREE(orbits, orbits_sz);

  nauty_free();

  return PermGroup(BSGS(nauty_generator_degree, nauty_generators, options));
}

} // namespace mpsym
