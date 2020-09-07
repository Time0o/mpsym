#ifndef GUARD_NAUTY_GRAPH_H
#define GUARD_NAUTY_GRAPH_H

#include <string>
#include <utility>
#include <vector>

#include "perm_group.hpp"
#include "perm_set.hpp"

namespace mpsym
{

namespace internal
{

class NautyGraph
{
public:
  NautyGraph(int n, int n_reduced, bool directed);

  ~NautyGraph();

  std::string to_gap() const;

  void add_edge(int from, int to);
  void set_partition(std::vector<std::vector<int>> const &ptn);

  PermGroup automorphisms(AutomorphismOptions const *options = nullptr) const;

private:
  unsigned long *_g;
  bool _directed;
  int _n, _n_reduced, _m;
  int *_lab, *_ptn, *_orbits;

  std::vector<std::pair<int, int>> _edges;
  std::vector<std::vector<int>> _ptn_expl;
};

} // namespace internal

} // namespace mpsym

#endif // GUARD_NAUTY_GRAPH_H
