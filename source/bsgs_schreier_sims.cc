#include <cassert>
#include <memory>
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
  DBG(DEBUG) << "Executing schreier sims algorithm";

  _strong_generators = generators;

  // intialize
  std::vector<PermSet> strong_generators;
  std::vector<std::vector<unsigned>> fundamental_orbits;
  std::vector<SchreierGeneratorQueue> schreier_generator_queues;

  schreier_sims_init(&strong_generators,
                     &fundamental_orbits,
                     &schreier_generator_queues);

  // main loop
  unsigned i = base_size();
  while (i >= 1u) {
top:
    schreier_generator_queues[i - 1].update(strong_generators[i - 1],
                                            fundamental_orbits[i - 1],
                                            _schreier_structures[i - 1]);

    for (Perm const &schreier_generator : schreier_generator_queues[i - 1]) {
      DBG(TRACE) << "Schreier Generator: " << schreier_generator;

      // strip
      TIMER_START("strip");

      auto strip_result(strip(schreier_generator, i));

      Perm strip_perm = std::get<0>(strip_result);
      unsigned strip_level = std::get<1>(strip_result);

      DBG(TRACE) << "Strips to: " << strip_perm << ", " << strip_level;

      TIMER_STOP("strip");

      if (strip_level < base_size() - i || !strip_perm.id()) {
        bool do_extend_base = i == base_size();

        if (do_extend_base) {
          TIMER_START("extend base");

          unsigned bp = 1u;
          for (;;) {
            auto it = std::find(_base.begin(), _base.end(), bp);

            if (it == _base.end() && strip_perm[bp] != bp)
              break;

            ++bp;

            assert(bp <= degree());
          }

          extend_base(bp);

          DBG(TRACE) << "Adjoined new basepoint:";
          DBG(TRACE) << "B = " << _base;

          strong_generators.emplace_back();
          fundamental_orbits.emplace_back();

          TIMER_STOP("extend base");
        }

        // update strong generators and fundamental orbits
        TIMER_START("update strong gens");

        DBG(TRACE) << "Updating strong generators:";

        // add result of 'strip' to strong generating set
        strong_generators[i].insert(strip_perm);

        // redetermine schreier structure and fundamental orbit
        update_schreier_structure(i, strong_generators[i]);
        fundamental_orbits[i] = orbit(i);

        DBG(TRACE) << "S(" << i + 1 << ") = " << strong_generators[i];
        DBG(TRACE) << "O(" << i + 1 << ") = " << fundamental_orbits[i];

        TIMER_STOP("update strong gens");

        // update schreier generator queue
        if (do_extend_base) {
          schreier_generator_queues.emplace_back();
        } else {
          schreier_generator_queues[i].invalidate();
        }

        ++i;
        goto top;
      }
    }
    --i;
  }

  schreier_sims_finish();
}

void BSGS::schreier_sims_random(PermSet const &generators, unsigned w)
{
  DBG(DEBUG) << "Executing (random) schreier sims algorithm";

  generators.assert_not_empty();

  _strong_generators = generators;

  // intialize
  std::vector<PermSet> strong_generators;
  std::vector<std::vector<unsigned>> fundamental_orbits;

  schreier_sims_init(&strong_generators, &fundamental_orbits);

  // random group element generator
  PrRandomizer pr(_strong_generators);

  unsigned c = 0u;
  while (c < w) {
    Perm rand_perm = pr.next();
    DBG(TRACE) << "Random group element: " << rand_perm;

    bool update_strong_generators = false;

    auto strip_result(strip(rand_perm));

    Perm strip_perm = std::get<0>(strip_result);
    unsigned strip_level = std::get<1>(strip_result);

    DBG(TRACE) << "Strips to: " << strip_perm << ", " << strip_level;

    if (strip_level <= base_size()) {
      update_strong_generators = true;

    } else if (!strip_perm.id()) {
      update_strong_generators = true;

      for (unsigned bp = 1u; bp <= degree(); ++bp) {
        if (strip_perm[bp] != bp) {
          extend_base(bp);
          strong_generators.emplace_back();

          DBG(TRACE) << "Adjoined new basepoint:";
          DBG(TRACE) << "B = " << _base;

          break;
        }
      }
    }

    if (update_strong_generators) {
      DBG(TRACE) << "Updating strong generators:";

      strong_generators.resize(strip_level);
      fundamental_orbits.resize(strip_level);

      for (unsigned i = 1u; i < strip_level; ++i) {
        strong_generators[i].insert(strip_perm);
        update_schreier_structure(i, strong_generators[i]);
        fundamental_orbits[i] = orbit(i);

        DBG(TRACE) << "S(" << (i + 1u) << ") = " << strong_generators[i];
        DBG(TRACE) << "O(" << (i + 1u) << ") = " << fundamental_orbits[i];
      }

      c = 0u;

    } else { ++c; }
  }

  schreier_sims_finish();
}

void BSGS::schreier_sims_init(
  std::vector<PermSet> *strong_generators,
  std::vector<std::vector<unsigned>> *fundamental_orbits,
  std::vector<SchreierGeneratorQueue> *schreier_generator_queues)
{
  DBG(DEBUG) << "=== Input";
  DBG(DEBUG) << "B = " << _base;
  DBG(DEBUG) << "S = " << _strong_generators;

  // add initial base points
  auto it = _strong_generators.begin();

  while (it != _strong_generators.end()) {
    Perm gen = *it;

    if (gen.id()) {
      it = _strong_generators.erase(it);

    } else {
      ++it;

      if (gen.stabilizes(_base.begin(), _base.end())) {
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

  strong_generators->resize(base_size());
  fundamental_orbits->resize(base_size());

  if (schreier_generator_queues)
    schreier_generator_queues->resize(base_size());

  // calculate initial strong generator sets
  for (unsigned i = 0u; i < base_size(); ++i) {
    for (Perm const &gen : _strong_generators) {
      if (gen.stabilizes(_base.begin(), _base.begin() + i))
        (*strong_generators)[i].insert(gen);
    }

    update_schreier_structure(i, (*strong_generators)[i]);
    (*fundamental_orbits)[i] = orbit(i);
  }

  DBG(DEBUG) << "=== Initial values";
  DBG(DEBUG) << "B = " << _base;
  for (unsigned i = 0u; i < base_size(); ++i) {
    DBG(DEBUG) << "S(" << (i + 1u) << ") = " << (*strong_generators)[i];
    DBG(DEBUG) << "O(" << (i + 1u) << ") = " << (*fundamental_orbits)[i];
  }
}

void BSGS::schreier_sims_finish()
{
  _strong_generators.clear();

  for (auto const &st : _schreier_structures) {
    auto stabilizers(st->labels());
    _strong_generators.insert(stabilizers.begin(), stabilizers.end());
  }

  _strong_generators.make_unique();

  DBG(DEBUG) << "=== Result";
  DBG(DEBUG) << "B = " << _base;
  DBG(DEBUG) << "SGS = " << _strong_generators;
}

} // namespace cgtl
