#include <cassert>
#include <climits>
#include <cmath>
#include <ctime>
#include <random>
#include <unordered_set>
#include <vector>

#include <boost/math/special_functions/prime.hpp>
#include <boost/multiprecision/miller_rabin.hpp>

#include "orbits.h"
#include "perm.h"
#include "pr_randomizer.h"
#include "util.h"

/**
 * @file pr_randomizer.cc
 * @brief Implements `PrRandomizer`.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

PrRandomizer::PrRandomizer(PermSet const &generators,
                           unsigned n_generators,
                           unsigned iterations)
: _gens_orig(generators)
{
  generators.assert_not_empty();

  _gens.insert(Perm(generators.degree()));

  if (generators.size() >= n_generators) {
    _gens.insert(generators.begin(), generators.end());
    n_generators = generators.size();
  } else {
    while (_gens.size() < n_generators) {
      unsigned missing = n_generators - _gens.size();
      if (missing > generators.size()) {
        _gens.insert(generators.begin(), generators.end());
      } else {
        _gens.insert(generators.begin(), generators.begin() + missing);
        break;
      }
    }
  }

  for (unsigned i = 0u; i < iterations; ++i)
    next();
}

Perm PrRandomizer::next()
{
  static auto re(util::random_engine());

  std::uniform_int_distribution<> randbool(0, 1);
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

bool PrRandomizer::test_symmetric(double epsilon)
{
  if (!test_altsym(epsilon))
    return false;

  return !generators_even();
}

bool PrRandomizer::test_alternating(double epsilon)
{
  if (!test_altsym(epsilon))
    return false;

  return generators_even();
}

bool PrRandomizer::test_altsym(double epsilon)
{
  assert(epsilon > 0.0);

  assert(_gens_orig.degree() >= 8u);

  // build prime number lookup table
  static std::unordered_set<unsigned> prime_lookup;

  if (prime_lookup.empty()) {
    //for (auto i = 0u; i <= boost::math::max_prime; ++i)
    for (auto i = 0u; i <= 1000u; ++i)
      prime_lookup.insert(boost::math::prime(i));
  }

  // check whether group is even transitive
  auto orbit(Orbit::generate(1, _gens_orig));

  if (orbit.size() != _gens_orig.degree())
    return false;

  // determine number of random elements to be tested
  double d = std::log(2.0) / std::log(static_cast<double>(_gens_orig.degree()));
  if (_gens_orig.degree() <= 16u)
    d *= 0.23;
  else
    d *= 0.39;

  double iterations_lower_bound = -std::log(epsilon) / d;
  assert(iterations_lower_bound < static_cast<double>(UINT_MAX));

  unsigned iterations =
    static_cast<unsigned>(std::ceil(iterations_lower_bound));

  // test whether random element contains p-cycle
  unsigned p_lower_bound = _gens_orig.degree() / 2u;
  unsigned p_upper_bound = _gens_orig.degree() - 2u;

  for (unsigned i = 0u; i < iterations; ++i) {
    for (auto const &cycle : next().cycles()) {
      auto cycle_len = cycle.size();

      bool is_prime = cycle_len <= boost::math::max_prime ?
        prime_lookup.find(cycle_len) != prime_lookup.end() :
        boost::multiprecision::miller_rabin_test(cycle_len, 25);

      if (!is_prime)
        continue;

      if (cycle_len > p_lower_bound && cycle_len < p_upper_bound)
        return true;
    }
  }

  return false;
}

bool PrRandomizer::generators_even()
{
  for (auto const &gen : _gens_orig) {
    if (!gen.even()) {
      return false;
    }
  }

  return true;
}

} // namespace cgtl
