#include <cassert>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

#include "hash.hpp"
#include "task_mapping.hpp"
#include "task_mapping_orbit.hpp"
#include "util.hpp"

namespace mpsym
{

void TMO::IterationState::advance()
{
  if (exhausted())
    return;

  auto current_copy(*current);
  _unprocessed.erase(current);

  _processed.insert(_hash(current_copy));

  for (auto const &gen : *_generators) {
    TaskMapping next(current_copy.permuted(gen));

    if (_processed.find(_hash(next)) == _processed.end())
      _unprocessed.insert(next);
  }

  current = _unprocessed.begin();
}

bool TMO::IterationState::exhausted() const
{ return _unprocessed.empty(); }

void TMO::IterationState::init_hash(TaskMapping const &root)
{
  auto support(_generators->support());

  std::set<unsigned> support_set(support.begin(), support.end());
  for (unsigned task : root)
    support_set.insert(task);

  unsigned n = support_set.size();
  unsigned k = root.size();

  uint64_t orbit_size_limit = 1;
  for (unsigned i = 0u; i < k; ++i) {
    orbit_size_limit *= n;
    if (orbit_size_limit >= std::numeric_limits<hash_type>::max()) {
      _hash = container_hash_truncated;
      return;
    }
  }

  unsigned i = 0u;
  for (unsigned task : support_set)
    _hash_support_map[task] = i++;

  _hash = [&](TaskMapping const &mapping){ return perfect_hash(mapping); };
}

TMO::IterationState::hash_type TMO::IterationState::perfect_hash(
  TaskMapping const &mapping) const
{
  hash_type h = 0u;

  hash_type factor = 1u;
  for (unsigned task : mapping) {
    auto it(_hash_support_map.find(task));
    assert(it != _hash_support_map.end());

    h += it->second * factor;
    factor *= _hash_support_map.size();
  }

  return h;
}

TMO::IterationState::hash_type TMO::IterationState::container_hash_truncated(
  TaskMapping const &mapping)
{ return util::container_hash(mapping.begin(), mapping.end()); }

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
