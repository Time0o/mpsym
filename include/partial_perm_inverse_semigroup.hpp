#if 0
#ifndef GUARD_PARTIAL_PERM_INVERSE_SEMIGROUP_H
#define GUARD_PARTIAL_PERM_INVERSE_SEMIGROUP_H

#include <cstddef>
#include <ostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "eemp.hpp"
#include "partial_perm.hpp"
#include "perm_group.hpp"
#include "util.hpp"

namespace mpsym
{

namespace internal
{

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
  PartialPermInverseSemigroup();
  PartialPermInverseSemigroup(std::vector<PartialPerm> const &generators);

  std::vector<PartialPerm> generators() const
  { return _generators; }

  void adjoin_generators(std::vector<PartialPerm> const &generators,
                         bool minimize = false);

  bool is_trivial() const
  { return _trivial; }

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

} // namespace internal

} // namespace mpsym

#endif // GUARD_PARTIAL_PERM_INVERSE_SEMIGROUP_H
#endif
