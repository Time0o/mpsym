#ifndef _GUARD_SCHREIER_TREE_H
#define _GUARD_SCHREIER_TREE_H

#include <map>
#include <vector>

#include "perm.h"

namespace cgtl
{

namespace schreier_sims
{

struct SchreierTree
{
  SchreierTree(unsigned degree) : _degree(degree) {}

  void create_root(unsigned root) { _root = root; }

  void create_edge(unsigned origin, unsigned destination, Perm const &perm) {
    _edges[origin] = destination;
    _labels[origin] = perm;
  }

  unsigned root() const { return _root; }

  std::vector<unsigned> nodes() const {
    std::vector<unsigned> result {_root};

    for (auto const &node : _edges)
      result.push_back(node.first);

    return result;
  }

  Perm transversal(unsigned origin) const
  {
    Perm result(_degree);

    unsigned current = origin;
    while(current != _root) {
      Perm const &label = _labels.find(current)->second;
      result = label * result;
      current = _edges.find(current)->second;
    }

    return result;
  }

  bool contains(unsigned node) const {
    return (node == _root) || (_edges.find(node) != _edges.end());
  }

private:
  unsigned _degree;
  unsigned _root = 0;
  std::map<unsigned, unsigned> _edges;
  std::map<unsigned, Perm> _labels;
};

} // namespace schreier_sims

} // namespace cgtl

#endif // _GUARD_SCHREIER_TREE_H
