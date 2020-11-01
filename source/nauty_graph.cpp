#include <cassert>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
  #include "nauty.h"
  #include "nausparse.h"
  #include "nautinv.h"
}

#include "arch_graph_system.hpp"
#include "bsgs.hpp"
#include "dump.hpp"
#include "nauty_graph.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"

namespace
{

template<typename T>
T *_alloc(std::size_t sz)
{
  void *ret = ALLOCS(sz, sizeof(T));
  if (!ret)
    throw std::bad_alloc();

  return static_cast<T *>(ret);
}

mpsym::internal::PermSet _gens;
int _gen_degree;

void _save_gens(int, int *perm, int *, int, int, int)
{
  std::vector<unsigned> tmp(_gen_degree);
  for (int i = 0; i < _gen_degree; ++i)
    tmp[i] = perm[i] + 1;

  _gens.emplace(tmp);
}

} // anonymous namespace

namespace mpsym
{

namespace internal
{

NautyGraph::NautyGraph(int n, int n_reduced, bool directed)
: _directed(directed),
  _n(n),
  _n_reduced(n_reduced)
{
#ifndef NDEBUG
  nauty_check(WORDSIZE, SETWORDSNEEDED(_n), _n, NAUTYVERSIONID);
#endif

  _lab = _alloc<int>(_n);
  _ptn = _alloc<int>(_n);
  _orbits = _alloc<int>(_n);
}

NautyGraph::~NautyGraph()
{
  FREES(_lab);
  FREES(_ptn);
  FREES(_orbits);

  naugraph_freedyn();
  nautil_freedyn();
  nauty_freedyn();
}

std::string NautyGraph::to_gap() const
{
  throw std::logic_error("TODO: need to consider loops");

#if 0
  std::stringstream ss;

  ss << "ReduceGroup(GraphAutoms([";

  // edges
  for (auto const &edge : _edges) {
    int source = edge.first + 1;
    int target = edge.second + 1;

    if (source != target)
      ss << "[" << source << "," << target << "],";
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
  ss << _n_reduced << ")";

  return ss.str();
#endif
}

void NautyGraph::add_edge(int from, int to)
{
  assert(from < _n);
  assert(to < _n);

  _edges.emplace_back(from, to);

  if (!_directed)
    _edges.emplace_back(to, from);
}

void NautyGraph::add_edges(std::map<int, std::vector<int>> const &adj)
{
  for (auto const &p : adj) {
    int from = p.first;
    for (int to : p.second)
      add_edge(from, to);
  }
}

void NautyGraph::set_partition(std::vector<std::vector<int>> const &ptn)
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

PermSet NautyGraph::automorphism_generators()
{
  if (_edges.empty())
    return {};

  // construct (sparse) nauty graph
  sparsegraph sg;

  SG_INIT(sg);

  int nde = static_cast<int>(_edges.size());

  SG_ALLOC(sg, _n, nde, "SG_ALLOC");

  sg.nv = _n;
  sg.nde = nde;

  std::vector<std::vector<int>> out_edges(_n);

  for (int v = 0; v < _n; ++v)
    sg.d[v] = 0;

  for (auto const &edge : _edges) {
    int source = edge.first;
    int target = edge.second;

    ++sg.d[source];

    out_edges[source].push_back(target);
  }

  int e_offs = 0u;
  for (int v = 0; v < _n; ++v) {
    sg.v[v] = e_offs;

    for (int target : out_edges[v])
      sg.e[e_offs++] = target;
  }

  // set nauty options
  static DEFAULTOPTIONS_SPARSEDIGRAPH(nauty_options_directed);
  static DEFAULTOPTIONS_SPARSEGRAPH(nauty_options_undirected);

  auto &nauty_options = _directed ? nauty_options_directed
                                  : nauty_options_undirected;

  nauty_options.defaultptn = _ptn_expl.empty() ? TRUE : FALSE;
  nauty_options.userautomproc = _save_gens;

  // call nauty
  _gens.clear();
  _gen_degree = _n_reduced;

  statsblk stats;
  sparsenauty(&sg, _lab, _ptn, _orbits, &nauty_options, &stats, nullptr);

  // free memory
  SG_FREE(sg);
  nausparse_freedyn();

  return _gens;
}

} // namespace internal

} // namespace mpsym
