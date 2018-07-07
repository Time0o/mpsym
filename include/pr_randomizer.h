#ifndef _GUARD_PR_RANDOMIZER_H
#define _GUARD_PR_RANDOMIZER_H

#include <random>
#include <vector>

#include "perm.h"

namespace cgtl
{

class PrRandomizer
{
public:
  PrRandomizer(std::vector<Perm> const &generators,
    unsigned n_generators = 10, unsigned iterations = 20);

  Perm next();

private:
  std::vector<Perm> _gens;
  unsigned _iterations;
};

} // namespace cgtl

#endif // _GUARD_PR_RANDOMIZER_H
