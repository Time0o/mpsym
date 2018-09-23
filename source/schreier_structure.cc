#include <vector>

#include "perm.h"
#include "schreier_structure.h"

namespace cgtl
{

void SchreierTree::create_root(unsigned root) { _root = root; }

void SchreierTree::create_labels(
  std::vector<Perm> const &labels) { _labels = labels; }

void SchreierTree::create_edge(
  unsigned origin, unsigned destination, unsigned label)
{

  _edges[origin] = destination;
  _edge_labels[origin] = label;
}

unsigned SchreierTree::root() const { return _root; }

std::vector<unsigned> SchreierTree::nodes() const
{
  std::vector<unsigned> result {_root};

  for (auto const &node : _edges)
    result.push_back(node.first);

  return result;
}

std::vector<Perm> SchreierTree::labels() const { return _labels; }

bool SchreierTree::contains(unsigned node) const
{
  return (node == _root) || (_edges.find(node) != _edges.end());
}

Perm SchreierTree::transversal(unsigned origin) const
{
  Perm result(_degree);

  unsigned current = origin;
  while(current != _root) {
    Perm const &label = _labels[_edge_labels.find(current)->second];
    result = label * result;
    current = _edges.find(current)->second;
  }

  return result;
}

} // namespace cgtl
