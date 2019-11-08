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
  // TODO: empty generator set

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
  // TODO: empty generator set

  std::vector<std::vector<unsigned>> res;

  std::vector<int> processed(generators.degree() + 1u, 0);
  unsigned num_processed = 0u;

  unsigned x = 1u;

  for (;;) {
    auto orbit(orbit_of(x, generators));

    res.emplace_back(orbit);

    if ((num_processed += orbit.size()) == generators.degree())
      break;

    for (unsigned y : orbit)
      processed[y] = 1;

    while (processed[x])
      ++x;
  }

  return res;
}

} // namespace cgtl
