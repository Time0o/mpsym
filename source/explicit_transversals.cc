#include <ostream>
#include <vector>

#include "perm.h"
#include "explicit_transversals.h"

namespace mpsym
{

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

void ExplicitTransversals::dump(std::ostream &os) const
{
  os << "explicit transversals:\n";

  for (auto const &tr : _orbit)
    os << tr.first << ": " << tr.second << "\n";
}

} // namespace mpsym
