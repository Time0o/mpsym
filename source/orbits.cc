#include <algorithm>
#include <cassert>
#include <memory>
#include <numeric>
#include <unordered_set>
#include <vector>

#include "orbits.h"
#include "perm.h"
#include "perm_set.h"
#include "schreier_structure.h"

namespace mpsym
{

Orbit Orbit::generate(unsigned x,
                      PermSet const &generators_,
                      std::shared_ptr<SchreierStructure> ss)
{
  Orbit orbit{x};

  if (generators_.trivial())
    return orbit;

  auto extend_ = [&](PermSet const &generators){
    orbit.extend(generators, {x}, {x}, ss);
  };

  if (ss) {
    generators_.assert_inverses();

    extend_(generators_);

  } else {
    PermSet generators(generators_);
    generators.insert_inverses();

    extend_(generators);
  }

  return orbit;
}

bool Orbit::generated_by(unsigned x, PermSet const &generators_) const
{
  if (generators_.trivial())
    return size() == 1u && (*this)[0] == x;

  PermSet generators(generators_);
  generators.insert_inverses();

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

    for (Perm const &gen : generators) {
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
                   PermSet const &generators_new,
                   std::shared_ptr<SchreierStructure> ss)
{
  auto generators(generators_old);
  generators.insert(generators_new.begin(), generators_new.end());

  if (ss) {
    generators.assert_inverses();

    for (Perm const &gen_new : generators_new)
      ss->add_label(gen_new);

  } else {
    generators.insert_inverses();
  }

  std::vector<unsigned> stack;
  std::unordered_set<unsigned> done(begin(), end());

  for (Perm const &gen_new : generators_new) {
    for (unsigned x : *this) {
      unsigned y = gen_new[x];

      if (done.find(y) != done.end())
        stack.push_back(y);
    }
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
      unsigned y = generators[i][x];

      if (done.find(y) == done.end()) {
        done.insert(y);
        stack.push_back(y);

        push_back(y);

        if (ss)
          ss->create_edge(y, x, i);
      }
    }
  }
}

OrbitPartition::OrbitPartition(unsigned degree)
: _partition_indices(degree, -1)
{}

OrbitPartition::OrbitPartition(unsigned degree,
                               std::vector<Orbit> const &partitions)
: _partitions(partitions),
  _partition_indices(degree, -1)
{ update_partition_indices(); }

OrbitPartition::OrbitPartition(unsigned degree,
                               std::vector<int> const &partition_indices)
: _partition_indices(partition_indices)
{
#ifndef NDEBUG
  assert(partition_indices.size() == degree);
#else
  (void)degree;
#endif

  update_partitions();
}

OrbitPartition::OrbitPartition(unsigned degree, PermSet const &generators)
: _partition_indices(degree, -1)
{
  if (generators.trivial())
    return;

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

  update_partition_indices();
}

std::vector<OrbitPartition> OrbitPartition::split(
  OrbitPartition const &split) const
{
  assert(split._partition_indices.size() == _partition_indices.size());

  std::vector<int> new_partition_indices(_partition_indices.size(), -1);
  std::vector<int> current_partitions_indices(split.num_partitions(), 0);

  std::vector<std::vector<int>> split_partition_indices(split.num_partitions());
  for (auto &partition_indices : split_partition_indices)
    partition_indices.resize(_partition_indices.size());

  for (unsigned x = 1u; x <= _partition_indices.size(); ++x) {
    int i = split.partition_index(x);
    int j = partition_index(x);

    if (new_partition_indices[j] == -1)
      new_partition_indices[j] = current_partitions_indices[i]++;

    split_partition_indices[i][x - 1u] = new_partition_indices[j];
  }

  std::vector<OrbitPartition> res;
  for (auto const &partition_indices : split_partition_indices)
    res.emplace_back(partition_indices.size(), partition_indices);

  return res;
}

void OrbitPartition::remove_from_partition(unsigned x)
{
  int i = _partition_indices[x - 1u];

  if (i == -1)
    return;

  _partitions[i].erase(
    std::find(_partitions[i].begin(), _partitions[i].end(), x));

  _partition_indices[x - 1u] = -1;
}

void OrbitPartition::change_partition(unsigned x, int i)
{
  if (i < 0) {
    remove_from_partition(x);
    return;
  }

  _partition_indices[x - 1u] = i;

  for (auto p_it = _partitions.begin(); p_it != _partitions.end(); ++p_it) {
    auto e_it = std::find(p_it->begin(), p_it->end(), x);

    if (e_it != p_it->end()) {
      p_it->erase(e_it);

      if (p_it->empty())
        _partitions.erase(p_it);

      break;
    }
  }

  add_to_partition(x, i);
}

void OrbitPartition::add_to_partition(unsigned x, int i)
{
  if (i >= static_cast<int>(_partitions.size()) - 1)
    _partitions.resize(i + 1);

  _partitions[i].push_back(x);
}

void OrbitPartition::update_partitions()
{
  for (unsigned x = 1u; x <= _partition_indices.size(); ++x) {
    int i = partition_index(x);

    if (i == -1)
      continue;

    add_to_partition(x, i);
  }
}

void OrbitPartition::update_partition_indices()
{
  for (int i = 0; i < static_cast<int>(_partitions.size()); ++i) {
    for (unsigned x : _partitions[i])
      _partition_indices[x - 1u] = i;
  }
}

} // namespace mpsym
