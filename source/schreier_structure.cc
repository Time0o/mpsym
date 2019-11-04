#include <vector>

#include "perm.h"
#include "schreier_structure.h"

namespace cgtl
{

void SchreierTree::create_root(unsigned root)
{
  _root = root;
}

void SchreierTree::create_labels(PermSet const &labels)
{
  _labels = labels;
}

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

PermSet SchreierTree::labels() const
{
  return _labels;
}

bool SchreierTree::contains(unsigned node) const
{
  return (node == _root) || (_edges.find(node) != _edges.end());
}

bool SchreierTree::incoming(unsigned node, Perm const &edge) const
{
  assert(edge.degree() == _degree);

  auto it = _edge_labels.find(edge[node]);
  if (it == _edge_labels.end())
    return false;

  return _labels[it->second] == edge;
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

void ExplicitTransversals::create_root(unsigned root)
{
  _root = root;
  _orbit[root] = Perm(_degree);
}

void ExplicitTransversals::create_labels(PermSet const &labels)
{
  _labels = labels;
}

void ExplicitTransversals::create_edge(
  unsigned origin, unsigned destination, unsigned label)
{
  if (_orbit.find(destination) == _orbit.end()) {
    _orbit[destination] = Perm(_degree);
    _orbit[origin] = _labels[label];
  } else {
    _orbit[origin] = _orbit[destination] * _labels[label];
  }
}

unsigned ExplicitTransversals::root() const
{
  return _root;
}

std::vector<unsigned> ExplicitTransversals::nodes() const
{
  std::vector<unsigned> res;
  for (auto item : _orbit)
    res.push_back(item.first);

  return res;
}

PermSet ExplicitTransversals::labels() const
{
  return _labels;
}

bool ExplicitTransversals::contains(unsigned node) const
{
  return _orbit.find(node) != _orbit.end();
}

bool ExplicitTransversals::incoming(unsigned, Perm const &) const
{
  return false;
}

Perm ExplicitTransversals::transversal(unsigned origin) const
{
  auto it(_orbit.find(origin));

  return it->second;
}

} // namespace cgtl
