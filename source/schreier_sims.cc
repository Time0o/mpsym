#include <cassert>
#include <memory>
#include <set>
#include <unordered_set>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "pr_randomizer.h"
#include "schreier_generator_queue.h"
#include "schreier_sims.h"
#include "schreier_structure.h"

/**
 * @file schreier_sims.cc
 * @brief Implements data structures and functions defined in schreier_sims.h.
 *
 * @author Timo Nicolai
 */

namespace
{

using cgtl::BSGS;
using cgtl::Perm;
using cgtl::schreier_sims::orbit;
using cgtl::schreier_sims::SchreierGeneratorQueue;

void schreier_sims_init(
  BSGS &bsgs,
  std::vector<std::vector<Perm>> &strong_generators,
  std::vector<std::vector<unsigned>> &fundamental_orbits,
  std::vector<SchreierGeneratorQueue> *schreier_generator_queues = nullptr)
{
  unsigned degree = bsgs.strong_generators[0].degree();

#ifndef NDEBUG
  for (auto const &gen : bsgs.strong_generators)
    assert(gen.degree() == degree && "all generators have same degree");
#endif

  Dbg(Dbg::DBG) << "=== Input";
  Dbg(Dbg::DBG) << "B = " << bsgs.base;
  Dbg(Dbg::DBG) << "S = " << bsgs.strong_generators;

  // add initial base points
  auto it  = bsgs.strong_generators.begin();

  while (it != bsgs.strong_generators.end()) {
    Perm gen = *it;

    if (gen.id()) {
      it = bsgs.strong_generators.erase(it);

    } else {
      ++it;

      bool stabilizes = true;
      for (unsigned b : bsgs.base) {
        if (gen[b] != b) {
          stabilizes = false;
          break;
        }
      }

      if (stabilizes) {
#ifndef NDEBUG
        bool extended_base = false;
#endif
        for (unsigned bp = 1u; bp <= degree; ++bp) {
          if (gen[bp] != bp) {
            bsgs.extend_base(bp);
#ifndef NDEBUG
            extended_base = true;
#endif
            break;
          }
        }

        assert(extended_base && "no generator fixes all base elements");
      }
    }
  }

  strong_generators.resize(bsgs.base.size());
  fundamental_orbits.resize(bsgs.base.size());

  // calculate initial strong generator sets
  for (unsigned i = 0u; i < bsgs.base.size(); ++i) {
    for (Perm const &gen : bsgs.strong_generators) {
      bool stabilizes = true;

      for (unsigned k = 0u; k < i; ++k) {
        if (gen[bsgs.base[k]] != bsgs.base[k]) {
          stabilizes = false;
          break;
        }
      }

      if (stabilizes)
        strong_generators[i].push_back(gen);
    }

    fundamental_orbits[i] = orbit(bsgs.base[i],
                                  strong_generators[i],
                                  bsgs.schreier_structures[i]);

    if (schreier_generator_queues)
      schreier_generator_queues->emplace_back(strong_generators[i],
                                              fundamental_orbits[i],
                                              bsgs.schreier_structures[i]);
  }

  Dbg(Dbg::DBG) << "=== Initial values";
  Dbg(Dbg::DBG) << "B = " << bsgs.base;
  for (unsigned i = 0u; i < bsgs.base.size(); ++i) {
    Dbg(Dbg::DBG) << "S(" << (i + 1u) << ") = " << strong_generators[i];
    Dbg(Dbg::DBG) << "O(" << (i + 1u) << ") = " << fundamental_orbits[i];
  }
}

void schreier_sims_finish(BSGS &bsgs)
{
  std::unordered_set<Perm> unique_generators;

  for (auto const &st : bsgs.schreier_structures) {
    auto stabilizers(st->labels());
    unique_generators.insert(stabilizers.begin(), stabilizers.end());
  }

  bsgs.strong_generators =
    std::vector<Perm>(unique_generators.begin(), unique_generators.end());

  orbit(bsgs.base[0], bsgs.strong_generators, bsgs.schreier_structures[0]);

  Dbg(Dbg::DBG) << "=== Result";
  Dbg(Dbg::DBG) << "B = " << bsgs.base;
  Dbg(Dbg::DBG) << "SGS = " << bsgs.strong_generators;
}

} // namespace

namespace cgtl
{

namespace schreier_sims
{

std::vector<unsigned> orbit(
  unsigned alpha, std::vector<Perm> const &generators,
  std::shared_ptr<SchreierStructure> st)
{
  assert(generators.size() > 0u && "generator set not empty");
  assert(alpha <= generators[0].degree() && "alpha <= N");

  std::vector<unsigned> result {alpha};

  if (st) {
    st->create_root(alpha);
    st->create_labels(generators);
  }

  std::vector<unsigned> stack {alpha};
  std::set<unsigned> done {alpha};

  while (!stack.empty()) {
    unsigned beta = stack.back();
    stack.pop_back();

    for (auto i = 0u; i < generators.size(); ++i) {
      Perm const &gen = generators[i];
      unsigned beta_prime = gen[beta];

      if (done.find(beta_prime) == done.end()) {
        result.push_back(beta_prime);
        done.insert(beta_prime);
        stack.push_back(beta_prime);

        if (st)
          st->create_edge(beta_prime, beta, i);
      }
    }
  }

  return result;
}

void schreier_sims(BSGS &bsgs)
{
  Dbg(Dbg::DBG) << "Executing schreier sims algorithm";

  assert(bsgs.strong_generators.size() > 0u && "generator set not empty");

  if (bsgs.strong_generators.size() == 1u && bsgs.strong_generators[0].id()) {
    bsgs.base.clear();
    return;
  }

  // intialize
  std::vector<std::vector<Perm>> strong_generators;
  std::vector<std::vector<unsigned>> fundamental_orbits;
  std::vector<SchreierGeneratorQueue> schreier_generator_queues;

  schreier_sims_init(bsgs,
                     strong_generators,
                     fundamental_orbits,
                     &schreier_generator_queues);

  // performance timers
  Timer_create("strip", Timer::MILLISECONDS);
  Timer_create("extend base", Timer::MILLISECONDS);
  Timer_create("update strong gens", Timer::MILLISECONDS);

  // main loop
  unsigned degree = bsgs.strong_generators[0].degree();

  unsigned i = bsgs.base.size();
  while (i >= 1u) {
top:
    schreier_generator_queues[i - 1].update(strong_generators[i - 1],
                                            fundamental_orbits[i - 1],
                                            bsgs.schreier_structures[i - 1]);

    for (Perm const &schreier_generator : schreier_generator_queues[i - 1]) {
      Dbg(Dbg::TRACE) << "Schreier Generator: " << schreier_generator;

      // strip
      Timer_start("strip");

      auto strip_result(bsgs.strip(schreier_generator, i));

      Perm strip_perm = std::get<0>(strip_result);
      unsigned strip_level = std::get<1>(strip_result);

      Dbg(Dbg::TRACE) << "Strips to: " << strip_perm << ", " << strip_level;

      Timer_stop("strip");

      if (strip_level < bsgs.base.size() - i || !strip_perm.id()) {
        bool extend_base = i == bsgs.base.size();

        if (extend_base) {
          Timer_start("extend base");

          unsigned bp = 1u;
          for (;;) {
            auto it = std::find(bsgs.base.begin(), bsgs.base.end(), bp);

            if (it == bsgs.base.end() && strip_perm[bp] != bp)
              break;

            assert(++bp <= degree);
          }

          bsgs.extend_base(bp);

          Dbg(Dbg::TRACE) << "Adjoined new basepoint:";
          Dbg(Dbg::TRACE) << "B = " << bsgs.base;

          strong_generators.emplace_back();
          fundamental_orbits.emplace_back();

          Timer_stop("extend base");
        }

        // update strong generators and fundamental orbits
        Timer_stop("update strong gens");

        Dbg(Dbg::TRACE) << "Updating strong generators:";

        // add result of 'strip' to strong generating set
        strong_generators[i].emplace_back(strip_perm);

        // redetermine schreier structure and fundamental orbit
        fundamental_orbits[i] = orbit(bsgs.base[i],
                                      strong_generators[i],
                                      bsgs.schreier_structures[i]);

        Dbg(Dbg::TRACE) << "S(" << i + 1 << ") = " << strong_generators[i];
        Dbg(Dbg::TRACE) << "O(" << i + 1 << ") = " << fundamental_orbits[i];

        Timer_stop("update strong gens");

        // update schreier generator queue
        if (extend_base) {
          schreier_generator_queues.emplace_back(strong_generators[i],
                                                 fundamental_orbits[i],
                                                 bsgs.schreier_structures[i]);
        } else {
          schreier_generator_queues[i].invalidate();
        }

        ++i;
        goto top;
      }
    }
    --i;
  }

  // dump performance statistics
  Timer_dump("strip");
  Timer_dump("extend base");
  Timer_dump("update strong gens");

  schreier_sims_finish(bsgs);
}

void schreier_sims_random(BSGS &bsgs, unsigned w)
{
  Dbg(Dbg::DBG) << "Executing (random) schreier sims algorithm";

  assert(bsgs.strong_generators.size() > 0u && "generator set not empty");

  if (bsgs.strong_generators.size() == 1u && bsgs.strong_generators[0].id()) {
    bsgs.base.clear();
    return;
  }

  // intialize
  std::vector<std::vector<Perm>> strong_generators;
  std::vector<std::vector<unsigned>> fundamental_orbits;

  schreier_sims_init(bsgs, strong_generators, fundamental_orbits);

  unsigned degree = bsgs.strong_generators[0].degree();

  // random group element generator
  PrRandomizer pr(bsgs.strong_generators);

  unsigned c = 0u;
  while (c < w) {
    Perm rand_perm = pr.next();
    Dbg(Dbg::TRACE) << "Random group element: " << rand_perm;

    bool update_strong_generators = false;

    auto strip_result(bsgs.strip(rand_perm));

    Perm strip_perm = std::get<0>(strip_result);
    unsigned strip_level = std::get<1>(strip_result);

    Dbg(Dbg::TRACE) << "Strips to: " << strip_perm << ", " << strip_level;

    if (strip_level <= bsgs.base.size()) {
      update_strong_generators = true;

    } else if (!strip_perm.id()) {
      update_strong_generators = true;

      for (unsigned bp = 1u; bp <= degree; ++bp) {
        if (strip_perm[bp] != bp) {
          bsgs.extend_base(bp);
          strong_generators.push_back({});

          Dbg(Dbg::TRACE) << "Adjoined new basepoint:";
          Dbg(Dbg::TRACE) << "B = " << bsgs.base;

          break;
        }
      }
    }

    if (update_strong_generators) {
      Dbg(Dbg::TRACE) << "Updating strong generators:";

      strong_generators.resize(strip_level);
      fundamental_orbits.resize(strip_level);

      for (unsigned i = 1u; i < strip_level; ++i) {
        strong_generators[i].push_back(strip_perm);
        fundamental_orbits[i] = orbit(bsgs.base[i], strong_generators[i],
                                      bsgs.schreier_structures[i]);

        Dbg(Dbg::TRACE) << "S(" << (i + 1u) << ") = " << strong_generators[i];
        Dbg(Dbg::TRACE) << "O(" << (i + 1u) << ") = " << fundamental_orbits[i];
      }

      c = 0u;

    } else { ++c; }
  }

  schreier_sims_finish(bsgs);
}

} // namespace schreier_sims

} // namespace cgtl
