#ifndef GUARD_PR_RANDOMIZER_H
#define GUARD_PR_RANDOMIZER_H

#include "perm_set.hpp"

namespace mpsym
{

namespace internal
{

class Perm;

class PrRandomizer
{
public:
  PrRandomizer(PermSet const &generators,
               unsigned n_generators = 10,
               unsigned iterations = 20);

  Perm next();

  bool test_symmetric(double epsilon = 1e-6);

private:
  bool test_altsym(double epsilon);
  bool generators_even();

  PermSet _gens_orig;
  PermSet _gens;
};

} // namespace internal

} // namespace mpsym

#endif // GUARD_PR_RANDOMIZER_H
