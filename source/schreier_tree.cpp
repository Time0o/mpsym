#include <ostream>
#include <utility>
#include <vector>

#include "perm.hpp"
#include "schreier_tree.hpp"

namespace mpsym
{

namespace internal
{

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

void SchreierTree::dump(std::ostream &os) const
{
  std::vector<std::vector<std::pair<unsigned, unsigned>>> adj(_degree + 1u);

  for (auto const &e : _edges) {
    auto label_it = _edge_labels.find(e.first);

    adj[e.first].emplace_back(e.second, label_it->second);
  }

  os << "schreier tree: [\n";

  for (auto origin = 1u; origin <= _degree; ++ origin) {
    if (adj[origin].empty())
      continue;

    os << "  " << origin << ": [";

    for (auto i = 0u; i < adj[origin].size(); ++i) {
      auto e = adj[origin][i];

      os << e.first << " " << _labels[e.second];

      if (i < adj[origin].size() - 1u)
        os << ", ";
    }

    os << "]\n";
  }

  os << "]\n";
}

} // namespace internal

} // namespace mpsym
