#ifndef _GUARD_PARTIAL_PERM_INVERSE_SEMIGROUP_H
#define _GUARD_PARTIAL_PERM_INVERSE_SEMIGROUP_H

#include <cstddef>
#include <ostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "eemp.h"
#include "partial_perm.h"
#include "perm_group.h"
#include "util.h"

/**
 * @file partial_perm_inverse_semigroup.h
 * @author Timo Nicolai
 *
 * @brief Defines `PartialPermInverseSemigroup`.
 */

namespace mpsym
{

/** A representation of inverse semigroups of partial permutations.
 *
 * Compared to PermGroup this class only supports very limited operations,
 * largely due to the lack of efficient algorithms to perform more complex
 * ones. This class nevertheless enables representation of inverse semigroups
 * of partial permutations without the need to explicitly store all elements,
 * allows for membership testing and extension of existing inverse semigroups
 * with additional generating elements.
 */
class PartialPermInverseSemigroup
{
  struct SccRepr
  {
    SccRepr() {}

    SccRepr(unsigned i,
            eemp::SchreierTree const &spanning_tree,
            PermGroup const &schreier_generators)
      : i(i),
        spanning_tree(spanning_tree),
        schreier_generators(schreier_generators)
      {}

    unsigned i;
    eemp::SchreierTree spanning_tree;
    PermGroup schreier_generators;
  };

public:
  /** Construct a trivial inverse semigroup of partial permutations.
   *
   * The trivial inverse semigroup of partial permutations contains as its only
   * element the unique partial permutation mapping from the empty set to itself
   */
  PartialPermInverseSemigroup();

  /** Construct an inverse semigroup of partial permutations from a set of
   *  generators.
   *
   * \param generators
   *     a generating set for the inverse semigroup of partial permutations to
   *     be constructed in the form of a vector of partial permutations, if this
   *     is empty, a trivial partial inverse semigroup of partial permutations
   *     is constructed
   */
  PartialPermInverseSemigroup(std::vector<PartialPerm> const &generators);

  /** Return a generating set for an inverse semigroup of partial permutations.
   *
   * Note that the returned generators are equal to those passed to this
   * object's constructor plus those added via adjoin() (unless the `minimize`
   * argument was set to `true`, in which case redundant generators passed via
   * adjoin() will not be included in this function's return value).
   *
   * \return a genering set of partial permutations in form of a vector
   *         containing the elements described above
   */
  std::vector<PartialPerm> generators() const
  { return _generators; }

  /** Adjoin additional partial permutations to the generating set of an inverse
   *  semigroup of partial permutations
   *
   * \param generators generators to be adjoined
   *
   * \param minimize
   *     if this is set to `true`, the generators are adjoined one by one and
   *     if in this process a generator is discovered to already be contained
   *     in the current inverse semigroup, it is skipped completely
   */
  void adjoin_generators(std::vector<PartialPerm> const &generators,
                         bool minimize = false);

  /** Check whether an inverse semigroup of partial permutations is trivial, i.e
   *  contains only the unique partial permutation mapping from the empty set to
   *  itself
   *
   * \return `true` if this inverse semigroup of partial permutations is
   *         trivial, else `false`
   */
  bool is_trivial() const
  { return _trivial; }

  /** Test membership of a partial permutation in an inverse semigroup of
   *  partial permutations.
   *
   * This function is guaranteed to be reasonably efficient, even for inverse
   * semigroups containing a very large number of elements.
   *
   * \param pperm
   *     the partial permutation to be tested for membership in this inverse
   *     semigroup
   *
   * \return `true` is `pperm` describes a partial permutation which is an
   *         element of the inverse semigroup described by this object
   */
  bool contains_element(PartialPerm const &pperm) const;

private:
  void update_action_component(std::vector<PartialPerm> const &generators);
  void update_scc_representatives();

  bool _trivial;

  std::vector<unsigned> _dom;
  std::vector<PartialPerm> _generators;

  std::vector<std::vector<unsigned>> _ac_im;
  std::unordered_map<std::vector<unsigned>,
                     unsigned,
                     util::ContainerHash<std::vector<unsigned>>> _ac_im_ht;
  eemp::SchreierTree _st_im;
  eemp::OrbitGraph _og_im;

  std::vector<unsigned> _scc;
  std::vector<SccRepr> _scc_repr;

  std::vector<PartialPerm> _r_class_repr;
};

std::ostream &operator<<(std::ostream &os,
                         PartialPermInverseSemigroup const &ppisg);

} // namespace mpsym

#endif // _GUARD_PARTIAL_PERM_INVERSE_SEMIGROUP_H
