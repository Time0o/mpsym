#ifndef _GUARD_SCHREIER_SIMS_H
#define _GUARD_SCHREIER_SIMS_H

#include <cassert>
#include <memory>
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

enum Construction {
  CONSTRUCTION_STANDARD,
  CONSTRUCTION_RANDOM,
  CONSTRUCTION_AUTO
};

enum Transversals {
  TRANSVERSALS_EXPLICIT,
  TRANSVERSALS_SCHREIER_TREES,
  TRANSVERSALS_SHALLOW_SCHREIER_TREES,
  TRANSVERSALS_AUTO
};

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
void schreier_sims(BSGS &bsgs);

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
void schreier_sims_random(BSGS &bsgs, unsigned w = 10);

} // namespace schreier_sims

} // namespace cgtl

#endif // _GUARD_SCHREIER_SIMS_H
