#include <cassert>
#include <stdexcept>
#include <utility>

#include "task_mapping.hpp"
#include "task_mapping_orbit.hpp"
#include "util.hpp"

namespace mpsym
{

void TMO::advance()
{
  if (exhausted())
    return;

  auto current_copy(*_current);
  _unprocessed.erase(_current);

  _processed.insert(current_copy);

  for (auto const &gen : _generators) {
    TaskMapping next(current_copy.permuted(gen));

    if (_processed.find(next) == _processed.end())
      _unprocessed.insert(next);
  }

  _current = _unprocessed.begin();
}

bool TMO::exhausted() const
{ return _unprocessed.empty(); }

std::pair<bool, unsigned> TMORs::insert(TaskMapping const &mapping)
{
  bool new_orbit;
  unsigned equivalence_class;

  auto it = _orbit_reprs.find(mapping);
  if (it == _orbit_reprs.end()) {
    new_orbit = true;
    equivalence_class = num_orbits();

    _orbit_reprs[mapping] = equivalence_class;

  } else {
    new_orbit = false;
    equivalence_class = it->second;
  }

  return {new_orbit, equivalence_class};
}

} // namespace mpsym
