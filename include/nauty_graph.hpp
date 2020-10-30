#ifndef GUARD_NAUTY_GRAPH_H
#define GUARD_NAUTY_GRAPH_H

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "perm_set.hpp"

namespace mpsym
{

namespace internal
{

class NautyGraph
{
public:
  NautyGraph(int n, bool directed)
  : NautyGraph(n, n, directed)
  {}

  NautyGraph(int n, int n_reduced, bool directed);

  ~NautyGraph();

  std::string to_gap() const;

  void add_edge(int from, int to);
  void add_edges(std::map<int, std::vector<int>> const &adj);

  void set_partition(std::vector<std::vector<int>> const &ptn);

  PermSet automorphism_generators();

private:
  bool _directed;
  int _n, _n_reduced;
  int *_lab, *_ptn, *_orbits;

  std::vector<std::pair<int, int>> _edges;
  std::vector<std::vector<int>> _ptn_expl;
};

} // namespace internal

} // namespace mpsym

#endif // GUARD_NAUTY_GRAPH_H
