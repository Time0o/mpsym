#ifndef _GUARD_ARCH_GRAPH_NAUTY_H
#define _GUARD_ARCH_GRAPH_NAUTY_H

#include <sstream>
#include <utility>
#include <vector>

#include "perm_set.h"

#include "dump.h"
#include "nauty.h"

namespace mpsym
{

namespace internal
{

class NautyGraph
{
public:
  NautyGraph(int n)
  : _n(n),
    _m(SETWORDSNEEDED(n))
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
    ADDONEEDGE(_g, from, to, _m);
    _edges.emplace_back(from, to);
  }

  void set_partition(std::vector<std::vector<int>> const &ptn)
  {
    __ptn = ptn;

    int i = 0;
    for (auto const &p : __ptn) {
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
        ss << "[" << source << "," << target << "],"
           << "[" << target << "," << source << "],";
      }
    }

    ss << "],";

    // partition
    auto ptn_inc(__ptn);
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
  int *_lab, *_ptn, *_orbits;

  std::vector<std::pair<int, int>> _edges;
  std::vector<std::vector<int>> __ptn;

  static PermSet _gens;
  static int _gen_degree;
};

} // namespace internal

} // namespace mpsym

#endif // _GUARD_ARCH_GRAPH_NAUTY_H
