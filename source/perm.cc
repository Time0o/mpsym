#include "perm.h"

#include <unordered_set>

namespace cgtl
{

std::unordered_set<unsigned> PermGroup::orbit(unsigned alpha) const
{
  assert(alpha <= _n);

  std::unordered_set<unsigned> result(_n);
  std::vector<unsigned> stack{alpha};

  while (!stack.empty()) {
    unsigned beta = stack.back();
    stack.pop_back();

    result.insert(beta);

    for (Perm const &gen : _generators) {
      unsigned beta_prime = gen[beta];

      if (result.find(beta_prime) == result.end())
        stack.push_back(beta_prime);
    }
  }

  return result;
}

}
