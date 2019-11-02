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
#include "schreier_structure.h"
#include "timer.h"

/**
 * @file bsgs_schreier_sims.cc
 * @brief Implements schreier sims algorithm variants.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

void BSGS::schreier_sims(PermSet const &generators)
{
  Dbg(Dbg::DBG) << "Executing schreier sims algorithm";

  strong_generators = generators.vect(); // TODO

  // intialize
  std::vector<PermSet> strong_generators_;
  std::vector<std::vector<unsigned>> fundamental_orbits;
  std::vector<SchreierGeneratorQueue> schreier_generator_queues;

  schreier_sims_init(&strong_generators_,
                     &fundamental_orbits,
                     &schreier_generator_queues);

  // performance timers
  Timer_create("strip", Timer::MILLISECONDS);
  Timer_create("extend base", Timer::MILLISECONDS);
  Timer_create("update strong gens", Timer::MILLISECONDS);

  // main loop
  unsigned i = base.size();
  while (i >= 1u) {
top:
    schreier_generator_queues[i - 1].update(strong_generators_[i - 1],
                                            fundamental_orbits[i - 1],
                                            schreier_structures[i - 1]);

    for (Perm const &schreier_generator : schreier_generator_queues[i - 1]) {
      Dbg(Dbg::TRACE) << "Schreier Generator: " << schreier_generator;

      // strip
      Timer_start("strip");

      auto strip_result(strip(schreier_generator, i));

      Perm strip_perm = std::get<0>(strip_result);
      unsigned strip_level = std::get<1>(strip_result);

      Dbg(Dbg::TRACE) << "Strips to: " << strip_perm << ", " << strip_level;

      Timer_stop("strip");

      if (strip_level < base.size() - i || !strip_perm.id()) {
        bool do_extend_base = i == base.size();

        if (do_extend_base) {
          Timer_start("extend base");

          unsigned bp = 1u;
          for (;;) {
            auto it = std::find(base.begin(), base.end(), bp);

            if (it == base.end() && strip_perm[bp] != bp)
              break;

            assert(++bp <= degree());
          }

          extend_base(bp);

          Dbg(Dbg::TRACE) << "Adjoined new basepoint:";
          Dbg(Dbg::TRACE) << "B = " << base;

          strong_generators_.emplace_back();
          fundamental_orbits.emplace_back();

          Timer_stop("extend base");
        }

        // update strong generators and fundamental orbits
        Timer_stop("update strong gens");

        Dbg(Dbg::TRACE) << "Updating strong generators:";

        // add result of 'strip' to strong generating set
        strong_generators_[i].insert(strip_perm);

        // redetermine schreier structure and fundamental orbit
        update_schreier_structure(i, strong_generators_[i]);
        fundamental_orbits[i] = orbit(i);

        Dbg(Dbg::TRACE) << "S(" << i + 1 << ") = " << strong_generators_[i];
        Dbg(Dbg::TRACE) << "O(" << i + 1 << ") = " << fundamental_orbits[i];

        Timer_stop("update strong gens");

        // update schreier generator queue
        if (do_extend_base) {
          schreier_generator_queues.emplace_back(strong_generators_[i],
                                                 fundamental_orbits[i],
                                                 schreier_structures[i]);
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

  schreier_sims_finish();
}

void BSGS::schreier_sims_random(PermSet const &generators, unsigned w)
{
  Dbg(Dbg::DBG) << "Executing (random) schreier sims algorithm";

  generators.assert_not_empty();

  strong_generators = generators.vect(); // TODO

  // intialize
  std::vector<PermSet> strong_generators_;
  std::vector<std::vector<unsigned>> fundamental_orbits;

  schreier_sims_init(&strong_generators_, &fundamental_orbits);

  // random group element generator
  PrRandomizer pr(strong_generators);

  unsigned c = 0u;
  while (c < w) {
    Perm rand_perm = pr.next();
    Dbg(Dbg::TRACE) << "Random group element: " << rand_perm;

    bool update_strong_generators = false;

    auto strip_result(strip(rand_perm));

    Perm strip_perm = std::get<0>(strip_result);
    unsigned strip_level = std::get<1>(strip_result);

    Dbg(Dbg::TRACE) << "Strips to: " << strip_perm << ", " << strip_level;

    if (strip_level <= base.size()) {
      update_strong_generators = true;

    } else if (!strip_perm.id()) {
      update_strong_generators = true;

      for (unsigned bp = 1u; bp <= degree(); ++bp) {
        if (strip_perm[bp] != bp) {
          extend_base(bp);
          strong_generators.emplace();

          Dbg(Dbg::TRACE) << "Adjoined new basepoint:";
          Dbg(Dbg::TRACE) << "B = " << base;

          break;
        }
      }
    }

    if (update_strong_generators) {
      Dbg(Dbg::TRACE) << "Updating strong generators:";

      strong_generators_.resize(strip_level);
      fundamental_orbits.resize(strip_level);

      for (unsigned i = 1u; i < strip_level; ++i) {
        strong_generators_[i].insert(strip_perm);
        update_schreier_structure(i, strong_generators_[i]);
        fundamental_orbits[i] = orbit(i);

        Dbg(Dbg::TRACE) << "S(" << (i + 1u) << ") = " << strong_generators_[i];
        Dbg(Dbg::TRACE) << "O(" << (i + 1u) << ") = " << fundamental_orbits[i];
      }

      c = 0u;

    } else { ++c; }
  }

  schreier_sims_finish();
}

void BSGS::schreier_sims_init(
  std::vector<PermSet> *strong_generators_,
  std::vector<std::vector<unsigned>> *fundamental_orbits,
  std::vector<SchreierGeneratorQueue> *schreier_generator_queues)
{
  Dbg(Dbg::DBG) << "=== Input";
  Dbg(Dbg::DBG) << "B = " << base;
  Dbg(Dbg::DBG) << "S = " << strong_generators;

  // add initial base points
  auto it = strong_generators.begin();

  while (it != strong_generators.end()) {
    Perm gen = *it;

    if (gen.id()) {
      it = strong_generators.erase(it);

    } else {
      ++it;

      bool stabilizes = true;
      for (unsigned b : base) {
        if (gen[b] != b) {
          stabilizes = false;
          break;
        }
      }

      if (stabilizes) {
#ifndef NDEBUG
        bool extended_base = false;
#endif
        for (unsigned bp = 1u; bp <= degree(); ++bp) {
          if (gen[bp] != bp) {
            extend_base(bp);
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

  strong_generators_->resize(base.size());
  fundamental_orbits->resize(base.size());

  // calculate initial strong generator sets
  for (unsigned i = 0u; i < base.size(); ++i) {
    for (Perm const &gen : strong_generators) {
      bool stabilizes = true;

      for (unsigned k = 0u; k < i; ++k) {
        if (gen[base[k]] != base[k]) {
          stabilizes = false;
          break;
        }
      }

      if (stabilizes)
        (*strong_generators_)[i].insert(gen);
    }

    update_schreier_structure(i, (*strong_generators_)[i]);
    (*fundamental_orbits)[i] = orbit(i);

    if (schreier_generator_queues)
      schreier_generator_queues->emplace_back((*strong_generators_)[i],
                                              (*fundamental_orbits)[i],
                                              schreier_structures[i]);
  }

  Dbg(Dbg::DBG) << "=== Initial values";
  Dbg(Dbg::DBG) << "B = " << base;
  for (unsigned i = 0u; i < base.size(); ++i) {
    Dbg(Dbg::DBG) << "S(" << (i + 1u) << ") = " << (*strong_generators_)[i];
    Dbg(Dbg::DBG) << "O(" << (i + 1u) << ") = " << fundamental_orbits[i];
  }
}

void BSGS::schreier_sims_finish()
{
  std::unordered_set<Perm> unique_generators;

  for (auto const &st : schreier_structures) {
    auto stabilizers(st->labels());
    unique_generators.insert(stabilizers.begin(), stabilizers.end());
  }

  strong_generators = PermSet(unique_generators.begin(),
                              unique_generators.end());

  update_schreier_structure(0, strong_generators);

  Dbg(Dbg::DBG) << "=== Result";
  Dbg(Dbg::DBG) << "B = " << base;
  Dbg(Dbg::DBG) << "SGS = " << strong_generators;
}

} // namespace cgtl
