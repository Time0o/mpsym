#include <cassert>
#include <ctime>
#include <random>
#include <vector>

#include "perm.h"
#include "pr_randomizer.h"

namespace cgtl
{

PrRandomizer::PrRandomizer(std::vector<Perm> const &generators,
  unsigned n_generators, unsigned iterations) : _iterations(iterations)
{
  assert(generators.size() > 0u &&
    "product replacement initialized with more than zero generators");

#ifndef NDEBUG
  for (auto const &gen : generators)
    assert(gen.degree() == generators[0].degree() &&
      "generators have same degree");
#endif

  _gens.push_back(Perm(generators[0].degree()));

  if (generators.size() >= n_generators) {
    _gens.insert(_gens.end(), generators.begin(), generators.end());
    n_generators = generators.size();
  } else {
    while (_gens.size() < n_generators) {
      unsigned missing = n_generators - _gens.size();
      if (missing > generators.size()) {
        _gens.insert(
          _gens.end(), generators.begin(), generators.end());
      } else {
        _gens.insert(
          _gens.end(), generators.begin(), generators.begin() + missing);
        break;
      }
    }
  }

  for (unsigned i = 0u; i < iterations; ++i)
    next();
}

Perm PrRandomizer::next()
{
  static std::default_random_engine re(time(nullptr));
  static std::uniform_int_distribution<> randbool(0, 1);

  std::uniform_int_distribution<> rands(1, _gens.size() - 1);
  std::uniform_int_distribution<> randt(1, _gens.size() - 1);

  int s, t;

  s  = rands(re);
  do { t = randt(re); } while (t == s);

  if (randbool(re)) {
    _gens[s] *= (randbool(re) ? _gens[t] : ~_gens[t]);
    _gens[0] *= _gens[s];
  } else {
    _gens[s] = (randbool(re) ? _gens[t] : ~_gens[t]) * _gens[s];
    _gens[0] = _gens[s] * _gens[0];
  }

  return _gens[0];
}

} // namespace cgtl
