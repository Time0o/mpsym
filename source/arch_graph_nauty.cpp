#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>

extern "C" {
  #include "nauty.h"
}

#include "arch_graph.hpp"
#include "arch_graph_system.hpp"
#include "bsgs.hpp"
#include "dump.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"

namespace mpsym
{

namespace internal
{

class NautyGraph
{
public:
  NautyGraph(int n, bool directed)
  : _n(n),
    _m(SETWORDSNEEDED(n)),
    _directed(directed)
  {
#ifndef NDEBUG
    nauty_check(WORDSIZE, _m, _n, NAUTYVERSIONID);
#endif

    _g = _alloc<graph>(_m * _n);
    _lab = _alloc<int>(_n);
    _ptn = _alloc<int>(_n);
    _orbits = _alloc<int>(_n);

    EMPTYGRAPH(_g, _m, _n);
  }

  ~NautyGraph()
  {
    FREES(_g);
    FREES(_lab);
    FREES(_ptn);
    FREES(_orbits);

    naugraph_freedyn();
    nautil_freedyn();
    nauty_freedyn();
  }

  void add_edge(int from, int to)
  {
    if (_directed)
      ADDONEARC(_g, from, to, _m);
    else
      ADDONEEDGE(_g, from, to, _m);

    _edges.emplace_back(from, to);
  }

  void set_partition(std::vector<std::vector<int>> const &ptn)
  {
    _ptn_expl = ptn;

    int i = 0;
    for (auto const &p : _ptn_expl) {
      for (auto j = 0u; j < p.size(); ++j) {
        _lab[i] = p[j];
        _ptn[i] = (j == p.size() - 1u) ? 0 : 1;
        ++i;
      }
    }
  }

  PermSet automorphisms(int reduce = 0) const
  {
    static DEFAULTOPTIONS_GRAPH(options);
    options.defaultptn = FALSE;
    options.userautomproc = _save_gens;

    _gens.clear();
    _gen_degree = reduce == 0 ? _n : reduce;

    statsblk stats;
    densenauty(_g, _lab, _ptn, _orbits, &options, &stats, _m, _n, nullptr);

    return _gens;
  }

  std::string to_gap(int reduce = 0) const
  {
    std::stringstream ss;

    ss << "ReduceGroup(GraphAutoms([";

    // edges
    for (auto const &edge : _edges) {
      int source = edge.first + 1;
      int target = edge.second + 1;

      if (source != target) {
        ss << "[" << source << "," << target << "],";

        if (!_directed)
           ss << "[" << target << "," << source << "],";
      }
    }

    ss << "],";

    // partition
    auto ptn_inc(_ptn_expl);
    for (auto &p : ptn_inc) {
      for (auto &v : p)
        ++v;
    }

    ss << DUMP(ptn_inc) << ",";

    // number of vertices
    ss << _n << "),";

    // number of vertices to reduce to
    ss << (reduce == 0 ? _n : reduce) << ")";

    return ss.str();
  }

private:
  template<typename T>
  static T *_alloc(std::size_t sz)
  {
    void *ret = ALLOCS(sz, sizeof(T));
    if (!ret)
      throw std::bad_alloc();

    return static_cast<T *>(ret);
  }

  static void _save_gens(int, int *perm, int *, int, int, int)
  {
    std::vector<unsigned> tmp(_gen_degree);
    for (int i = 0; i < _gen_degree; ++i)
      tmp[i] = perm[i] + 1;

    _gens.emplace(tmp);
  }

  graph * _g;
  std::size_t _n, _m;
  bool _directed;
  int *_lab, *_ptn, *_orbits;

  std::vector<std::pair<int, int>> _edges;
  std::vector<std::vector<int>> _ptn_expl;

  static PermSet _gens;
  static int _gen_degree;
};

PermSet NautyGraph::_gens;
int NautyGraph::_gen_degree;

} // namespace internal

internal::NautyGraph ArchGraph::graph_nauty() const
{
  int cts = _channel_types.size();
  int cts_log2 = 0; while (cts >>= 1) ++cts_log2;

  int n_orig = num_processors();
  int n = n_orig * (cts_log2 + 1u);

  internal::NautyGraph g(n, _directed);

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

std::string ArchGraph::to_gap_nauty() const
{
  auto g(graph_nauty());

  return g.to_gap(num_processors());
}

internal::PermGroup ArchGraph::automorphisms_nauty(
  AutomorphismOptions const *options) const
{
  auto g(graph_nauty());

  auto automs(g.automorphisms(num_processors()));

  return internal::PermGroup(internal::BSGS(num_processors(), automs, options));
}

} // namespace mpsym
