#include <algorithm>
#include <cassert>
#include <memory>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "perm.hpp"
#include "perm_set.hpp"

namespace mpsym
{

namespace internal
{

unsigned PermSet::smallest_moved_point() const
{
  assert(!trivial());

  for (unsigned smp = 0u; smp < degree(); ++smp) {
    for (auto const &perm : *this) {
      if (perm[smp] != smp)
        return smp;
    }
  }

  throw std::logic_error("unreachable");
}

unsigned PermSet::largest_moved_point() const
{
  assert(!trivial());

  for (int lmp = static_cast<int>(degree() - 1u); lmp >= 0; --lmp) {
    unsigned lmp_ = static_cast<unsigned>(lmp);

    for (auto const &perm : *this) {
      if (perm[lmp_] != lmp_)
        return lmp_;
    }
  }

  throw std::logic_error("unreachable");
}

std::vector<unsigned> PermSet::support() const
{
  std::vector<unsigned> sup;

  for (unsigned x = 0u; x < degree(); ++x) {
    bool moved = false;
    for (auto const &perm : *this) {
      if (perm[x] != x) {
        moved = true;
        break;
      }
    }

    if (moved)
      sup.push_back(x);
  }

  return sup;
}

void PermSet::make_unique()
{
  std::vector<Perm> unique_perms;

  std::unordered_set<Perm> seen;
  for (Perm const &perm : _perms) {
    if (seen.find(perm) != seen.end())
      continue;

    unique_perms.push_back(perm);
    seen.insert(perm);
  }

  _perms = unique_perms;
}

void PermSet::insert_inverses()
{
  auto perms_and_inverses(_perms);

  for (auto const &perm : *this)
    perms_and_inverses.emplace_back(~perm);

  _perms = perms_and_inverses;

  make_unique();
}

void PermSet::minimize_degree()
{
  if (empty())
    return;

  std::vector<std::vector<unsigned>> moved_sets(size());

  std::vector<unsigned> compression_mapping(degree());
  for (unsigned i = 0u; i < degree(); ++i)
    compression_mapping[i] = i;

  std::queue<unsigned> non_moved_queue;

  unsigned new_degree = 1u;

  for (unsigned i = 0u; i < degree(); ++i) {
    bool moved = false;
    for (auto j = 0u; j < _perms.size(); ++j) {
      if (_perms[j][i] != i) {
        moved_sets[j].push_back(i);
        moved = true;
      }
    }

    if (moved) {
      if (!non_moved_queue.empty()) {
        unsigned compress_to = non_moved_queue.front();
        compression_mapping[i] = compress_to;
        new_degree = compress_to + 1u;

        non_moved_queue.pop();
        non_moved_queue.push(i);
      } else {
        new_degree = i + 1u;
      }
    } else {
      non_moved_queue.push(i);
    }
  }

  std::vector<unsigned> id(new_degree);
  std::iota(id.begin(), id.end(), 0u);

  for (unsigned i = 0u; i < _perms.size(); ++i) {
    auto gen(id);
    for (unsigned j = 0u; j < moved_sets[i].size(); ++j) {
      unsigned x = moved_sets[i][j];
      unsigned y = _perms[i][x];

      gen[compression_mapping[x]] = compression_mapping[y];
    }

    _perms[i] = Perm(gen);
  }
}

} // namespace internal

} // namespace mpsym
