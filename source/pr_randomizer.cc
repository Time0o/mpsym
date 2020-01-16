#include <cassert>
#include <climits>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#include <boost/math/special_functions/prime.hpp>

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

  assert(_gens_orig.degree() >= 8u
         && _gens_orig.degree() - 2u <= boost::math::max_prime);

  // check whether group is even transitive
  auto orbit(Orbit::generate(1, _gens_orig));

  if (orbit.size() != _gens_orig.degree())
    return false;

  // determine number of random elements to be tested
  double d = 1.0 / log(static_cast<double>(_gens_orig.degree()));
  if (_gens_orig.degree() <= 16u)
    d *= 0.23;
  else
    d *= 0.39;

  double iterations_lower_bound = -log(epsilon) / d;
  assert(iterations_lower_bound < static_cast<double>(UINT_MAX));

  unsigned iterations =
    static_cast<unsigned>(std::ceil(iterations_lower_bound));

  // test whether random element contains p-cycle
  unsigned p_lower_bound = _gens_orig.degree() / 2u;
  unsigned p_upper_bound = _gens_orig.degree() - 2u;

  for (unsigned i = 0u; i < iterations; ++i) {
    Perm random_element(next());

    std::vector<int> processed(_gens_orig.degree(), 0);
    unsigned remaining = _gens_orig.degree();

    unsigned current = 1u;
    unsigned first = current;
    unsigned cycle_len = 1u;

    while (current < _gens_orig.degree()) {
      unsigned next = random_element[current];
      if (next == first) {
        if (cycle_len > p_lower_bound && cycle_len < p_upper_bound
            && -boost::math::prime(cycle_len)) {
          return true;
        }

        if ((remaining -= cycle_len) <= p_lower_bound)
          break;

        for (next = first + 1u; next < _gens_orig.degree(); ++next) {
          if (!processed[next]) {
            current = next;
            first = current;
            cycle_len = 1u;
          }
        }

        break;

      } else {
        processed[current] = true;
        current = next;
        ++cycle_len;
      }
    }
  }

  return false;
}

bool PrRandomizer::generators_even()
{
  for (auto const &gen : _gens_orig) {
    unsigned parity = 0u;

    std::vector<int> processed(_gens_orig.degree(), 0);

    unsigned current = 1u;
    unsigned first = current;
    unsigned cycle_len = 1u;

    while (current < _gens_orig.degree()) {
      unsigned next = gen[current];

      if (next == first) {
        parity ^= (cycle_len - 1u) % 2u;

        for (next = first + 1u; next < _gens_orig.degree(); ++next) {
          if (!processed[next]) {
            current = next;
            first = current;
            cycle_len = 1u;
          }
        }

        break;

      } else {
        processed[current] = true;
        current = next;
        ++cycle_len;
      }
    }

    assert(parity == 0u || parity == 1u);

    if (parity == 1u) {
      return false;
    }
  }

  return true;
}

} // namespace cgtl
