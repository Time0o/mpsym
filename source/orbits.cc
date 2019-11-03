#include <set>
#include <vector>

#include "orbits.h"
#include "perm.h"
#include "perm_set.h"
#include "schreier_structure.h"

namespace cgtl
{

std::vector<unsigned>
orbit_of(unsigned x, PermSet const &generators, SchreierStructure *ss)
{
  std::vector<unsigned> res {x};

  if (ss) {
    ss->create_root(x);
    ss->create_labels(generators);
  }

  std::vector<unsigned> stack {x};
  std::set<unsigned> done {x};

  while (!stack.empty()) {
    unsigned beta = stack.back();
    stack.pop_back();

    for (auto i = 0u; i < generators.size(); ++i) {
      Perm gen = generators[i];
      unsigned beta_prime = gen[beta];

      if (done.find(beta_prime) == done.end()) {
        done.insert(beta_prime);
        stack.push_back(beta_prime);
        res.push_back(beta_prime);

        if (ss)
          ss->create_edge(beta_prime, beta, i);
      }
    }
  }

  return res;
}

std::vector<std::vector<unsigned>>
orbit_partition(PermSet const &generators)
{
  std::vector<std::vector<unsigned>> res;

  std::vector<int> orbit_indices(generators.degree() + 1u, -1);

  unsigned processed = 0u;

  for (auto i = 1u; i <= generators.degree(); ++i) {
    int orbit_index1 = orbit_indices[i];
    if (orbit_index1 == -1) {
      orbit_index1 = static_cast<int>(res.size());
      orbit_indices[i] = orbit_index1;

      res.push_back({i});

      if (++processed == generators.degree())
        return res;
    }

    for (Perm const &gen : generators) {
      unsigned j = gen[i];

      int orbit_index2 = orbit_indices[j];
      if (orbit_index2 == -1) {
        res[orbit_index1].push_back(j);
        orbit_indices[j] = orbit_index1;

        if (++processed == generators.degree())
          return res;
      }
    }
  }

  return res;
}

} // namespace cgtl
