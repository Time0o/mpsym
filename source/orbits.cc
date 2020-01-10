#include <algorithm>
#include <memory>
#include <unordered_set>
#include <vector>

#include "orbits.h"
#include "perm.h"
#include "perm_set.h"
#include "schreier_structure.h"

namespace cgtl
{

Orbit Orbit::generate(unsigned x,
                      PermSet const &generators,
                      std::shared_ptr<SchreierStructure> ss)
{
  Orbit orbit{x};

  if (!generators.empty())
    orbit.extend(generators, {x}, {x}, ss);

  return orbit;
}

bool Orbit::generated_by(unsigned x, PermSet const &generators) const
{
  if (generators.empty())
    return size() == 1u && (*this)[0] == x;

  std::vector<int> in_orbit_ref(generators.degree() + 1u, 0);
  std::vector<int> in_orbit(generators.degree() + 1u, 0);

  for (unsigned y : *this)
    in_orbit_ref[y] = 1;

  if (!in_orbit_ref[x])
    return false;

  in_orbit[x] = 1;

  std::vector<unsigned> stack{x};

  auto target = size() - 1u;

  while (!stack.empty()) {
    unsigned x = stack.back();
    stack.pop_back();

    for (auto const &gen : generators) {
      unsigned y = gen[x];

      if (!in_orbit_ref[y])
        return false;

      if (in_orbit[y] == 0) {
        in_orbit[y] = 1;

        stack.push_back(y);

        if (--target == 0u) {
          return true;
        }
      }
    }
  }

  return false;
}

void Orbit::update(PermSet const &generators_old,
                   Perm const &generator_new,
                   std::shared_ptr<SchreierStructure> ss)
{
  auto generators(generators_old);
  generators.insert(generator_new);

  if (ss)
    ss->add_label(generator_new);

  std::vector<unsigned> stack;
  std::unordered_set<unsigned> done(begin(), end());

  for (unsigned x : *this) {
    unsigned x_prime = generator_new[x];
    if (done.find(x_prime) != done.end())
      stack.push_back(x_prime);
  }

  extend(generators, stack, done, ss);
}

void Orbit::extend(PermSet const &generators,
                   std::vector<unsigned> stack,
                   std::unordered_set<unsigned> done,
                   std::shared_ptr<SchreierStructure> ss)
{
  while (!stack.empty()) {
    unsigned x = stack.back();
    stack.pop_back();

    for (auto i = 0u; i < generators.size(); ++i) {
      unsigned x_prime = generators[i][x];

      if (done.find(x_prime) == done.end()) {
        done.insert(x_prime);
        stack.push_back(x_prime);

        push_back(x_prime);

        if (ss)
          ss->create_edge(x_prime, x, i);
      }
    }
  }
}

OrbitPartition::OrbitPartition(PermSet const &generators)
: _partition_indices(generators.degree())
{
  generators.assert_not_empty();

  // determine partitions
  std::vector<int> processed(generators.degree() + 1u, 0);
  unsigned num_processed = 0u;

  unsigned x = 1u;

  for (;;) {
    auto orbit(Orbit::generate(x, generators));

    _partitions.push_back(orbit);

    if ((num_processed += orbit.size()) == generators.degree())
      break;

    for (unsigned y : orbit)
      processed[y] = 1;

    while (processed[x])
      ++x;
  }

  // determine partition indices
  for (auto i = 0u; i < _partitions.size(); ++i) {
    for (unsigned x : _partitions[i])
      _partition_indices[x - 1u] = i;
  }
}

} // namespace cgtl
