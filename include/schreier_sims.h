#ifndef _GUARD_SCHREIER_SIMS_H
#define _GUARD_SCHREIER_SIMS_H

#include <cassert>
#include <memory>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "pr_randomizer.h"
#include "schreier_generator_queue.h"
#include "schreier_sims.h"
#include "schreier_structure.h"
#include "timer.h"

/**
 * @file schreier_sims.h
 * @brief Definitions related to the Schreier-Sim algorithm.
 *
 * For details on the Schreier-Sims algorithm consult \cite holt05.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

namespace schreier_sims
{

/** Compute the orbit decomposition of a permutation group given as a set of
 *  generating permutations.
 *
 * The *orbit* of a group element \f$g \in G\f$ where \f$G\f$ acts on the set
 * \f$\Omega\f$ is the set \f$G(x) = \{x^g \mid g \in G\}\f$. This function
 * returns the set \f$\{\{x^g \mid g \in G\} \mid x \in \Omega\}\f$ which
 * contains every possible element orbit and is always a partition of \f$G\f$.
 *
 * \param generators
 *     a vector of `Perm` objects describing a generating set for the
 *     permutation group \f$G\f$
 *
 * \return \f$G\f$'s orbit partition as a vector of vectors where each element
 *         vector contains all elements in one of the possible element orbits
 *         (sorted in ascending order) and the element vectors themselves are
 *         sorted by their the smallest element (in ascending order)
 */
std::vector<std::vector<unsigned>> orbits(std::vector<Perm> const &generators);

/** Compute the orbit of a single element \f$x \in \Omega\f$ for the group \f$G
 *  \leq Sym(\Omega)\f$ described by a set of generating permutations.
 *
 * The *orbit* is the set \f$G(x) = \{x^g \mid g \in G\}\f$. This function also
 * computes the corresponding *orbit transversals*, i.e. the set \f$\{u_{\beta}
 * \mid \beta \in G(x)\}\f$ where \f$u_{\beta} \in G\f$ maps \f$x\f$ to
 * \f$\beta\f$. These orbit transversals are stored in a *Schreier structure*
 * (see schreier_sims() and SchreierStructure).
 *
 * \param x
 *     the element \f$x \in \Omega\f$ for which \f$G(x)\f$ should be computed
 *
 * \param generators
 *     a vector of `Perm` objects describing a generating set for the
 *     permutation group \f$G\f$
 *
 * \param[out] st
 *     a pointer to a preallocated Schreier structure which will store the orbit
 *     transversals after the execution of this function, provided
 *     `st != nullptr`
 *
 * \return the orbit \f$G(x)\f$ as described above in the form of an ordered
 *         vector of positive integers \f$\in \Omega\f$
 */
std::vector<unsigned> orbit(
  unsigned x, std::vector<Perm> const &generators,
  std::shared_ptr<SchreierStructure> st = nullptr);

/** Test membership of an arbitrary element in \f$Sym(\Omega)\f$ in a
 *  permutation group \f$G\f$ acting on \f$\Omega\f$.
 *
 * This function computes a pair containing a permutation \f$h\f$ and a
 * positive integer \f$i\f$. If the permutation \f$g\f$ described by `perm` is
 * contained in the group \f$G\f$ with *base* `base` and *base orbits/base
 * orbit transversals* described by the *Schreier structures* `sts` (see also
 * schreier_sims()), \f$h\f$ will be the identity permutation and \f$i\f$ will
 * be \f$k + 1\f$ (where \f$k\f$ is the length of the given base). This implies
 * that \f$g\f$ can be decomposed as \f$u_k u_{k-1} \dots u_i\f$ with \f$u_i
 * \in U^{(i)}\f$ where \f$U^{(i)}\f$ is the \f$i\f$th *base orbit transversal*.
 *
 * If \f$g \notin G\f$, this will not hold. The returned values \f$h\f$ and
 * \f$i\f$ are in this case nevertheless vital to the implementation of the
 * *Schreier-Sims* algorithm in schreier_sims() (This is an implementation
 * detail, refer to \cite holt05 for further information).
 *
 * \param perm
 *     the permutation which should be tested for membership in the group with
 *     base `base` and base orbits/base orbit transversals described by `sts`
 *
 * \param bsgs
 *     a base and strong generating set, e.g. constructed using schreier_sims()
 *
 * \param offs
 *     only consider base elements \f$\beta_i\f$ with \f$i >\f$ `offs`
 *
 * \return a pair containing a permutation \f$h\f$ and a positive integer
 *         \f$i\f$ as described above
 */
std::pair<Perm, unsigned> strip(Perm const &perm,
                                BSGS const &bsgs,
                                unsigned offs = 0);

template <typename SS>
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
        for (unsigned i = 1u; i <= degree; ++i) {
          if (gen[i] != i) {
            bsgs.base.push_back(i);
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

  for (auto i = bsgs.schreier_structures.size(); i < bsgs.base.size(); ++i)
    bsgs.schreier_structures.push_back(std::make_shared<SS>(degree));

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

void schreier_sims_finish(BSGS &bsgs);

/** The simple *Schreier-Sims* algorithm
 *
 * The following explanation is based on \cite holt05, consult the book for
 * a more detailed depiction.
 *
 * This algorithm computes a *base and strong generating set* (short *BSGS*)
 * completely representing a specific permutation group.  A base of a
 * permutation group \f$\f$ is a list \f$B = [\beta_1, \dots, \beta_k]\f$ of
 * elements \f$\beta_i \in \Omega,\ 1 \leq i \leq k\f$ (where \f$\Omega\f$ is
 * the set on which \f$G\f$ acts) where the only element in \f$G\f$ that fixes
 * every element in \f$B\f$ is the identity.
 *
 * A strong generating set is a set of partial permutations such that for \f$i
 * = 1, 2, \dots, k\f$ it holds that \f$G^{(i)} = \left<S^{(i)}\right>\f$.
 * Here, \f$G^{(i)}\f$ is the \f$i\f$th *base stabilizer* which maps all
 * elements \f$\beta_1, \dots, \beta_{i-1}\f$ to themselves and \f$S^{(i)} := S
 * \cap G^{(i)}\f$.
 *
 * This function also calculates *base orbits/base orbit transversals* in the
 * form of a vector of *Schreier structures* (one per base element). The
 * \f$i\f$th base orbit is the set \f$\Delta^{(i)} := \beta_i^{H^{(i)}}\f$
 * (where \f$H^{(i)} := \left<S^{(i)}\right>\f$) and the corresponding base
 * orbit transversal is \f$U^{(i)} := \{u_{\beta}^{(i)} \mid \beta \in
 * \Delta^{(i)}\}\f$ (where \f$u_{\beta}^{(i)} \in G\f$ maps \f$\beta_i\f$ to
 * \f$\beta\f$).
 *
 * Note that in general there is no unique BSGS for a given permutation group
 * and this function makes no guarantees as to which possible BSGS is returned.
 *
 * \param[in,out] bsgs
 *     initial base and generators in the form of a base and strong generating
 *     set (BSGS) object; the resulting base will be an extension of the
 *     initial base, i.e. in general more elements will be appended to it; the
 *     base may be empty initially (and should be if no particular base prefix
 *     is desired); the generators are replaced by the computed strong
 *     generating set (in no particular order); this parameter's Schreier
 *     structures, which describe the orbits of base elements and the associated
 *     transversals, are also set by this function
 */
template <typename SS>
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

  schreier_sims_init<SS>(bsgs,
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

      std::pair<Perm, unsigned> strip_result =
        strip(schreier_generator, bsgs, i);

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

          bsgs.base.emplace_back(bp);
          bsgs.schreier_structures.emplace_back(std::make_shared<SS>(degree));

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

/** The *Random Schreier-Sims* algorithm
 *
 * This is a variant of the algorithm implemented by schreier_sims() which uses
 * randomly generated group elements (via PrRandomizer) to obtain a BSGS. The
 * generated BSGS is not guaranteed to be correct BSGS for the permutation group
 * described by the given generators, to achieve certainty, an additional call
 * to TODO is necessary.
 *
 * \param[in,out] bsgs
 *     initial base and generators in the form of a base and strong generating
 *     set (BSGS) object; the resulting base will be an extension of the
 *     initial base, i.e. in general more elements will be appended to it; the
 *     base may be empty initially (and should be if no particular base prefix
 *     is desired); the generators are replaced by the computed strong
 *     generating set (in no particular order); this parameter's Schreier
 *     structures, which describe the orbits of base elements and the associated
 *     transversals, are also set by this function
 *
 * \param w
 *     tuning parameter, higher values increase the probability that the
 *     computed BSGS is actually a BSGS for the permutation group generated by
 *     the given generators; more specifically this probability is approximately
 *     \f$1 - 2^{-w}\f$, refer to \cite holt05 for details
 */
template <typename SS>
void schreier_sims_random(BSGS &bsgs, unsigned w = 10)
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

  schreier_sims_init<SS>(bsgs, strong_generators, fundamental_orbits);

  unsigned degree = bsgs.strong_generators[0].degree();

  // random group element generator
  PrRandomizer pr(bsgs.strong_generators);

  unsigned c = 0u;
  while (c < w) {
    Perm rand_perm = pr.next();
    Dbg(Dbg::TRACE) << "Random group element: " << rand_perm;

    bool update_strong_generators = false;
    std::pair<Perm, unsigned> strip_result = strip(rand_perm, bsgs);

    Perm strip_perm = std::get<0>(strip_result);
    unsigned strip_level = std::get<1>(strip_result);
    Dbg(Dbg::TRACE) << "Strips to: " << strip_perm << ", " << strip_level;

    if (strip_level <= bsgs.base.size()) {
      update_strong_generators = true;

    } else if (!strip_perm.id()) {
      update_strong_generators = true;

      for (unsigned i = 1u; i <= degree; ++i) {
        if (strip_perm[i] != i) {
          bsgs.base.push_back(i);
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
      bsgs.schreier_structures.resize(strip_level, std::make_shared<SS>(degree));

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

#endif // _GUARD_SCHREIER_SIMS_H
