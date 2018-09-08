#ifndef _GUARD_PARTIAL_PERM_INVERSE_SEMIGROUP_H
#define _GUARD_PARTIAL_PERM_INVERSE_SEMIGROUP_H

#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "boost/container_hash/hash.hpp"

#include "eemp.h"
#include "partial_perm.h"
#include "perm_group.h"

namespace cgtl
{

class PartialPermInverseSemigroup
{
  struct VectorHash {
    std::size_t operator()(std::vector<unsigned> const &v) const {
      return boost::hash_range(v.begin(), v.end());
    }
  };

  struct SccRepr {
    SccRepr() {}
    SccRepr(unsigned i,
            EEMP::OrbitGraph const &orbit_graph,
            std::vector<unsigned> const &scc,
            PermGroup const &schreier_generators);

    unsigned i;
    std::map<unsigned, unsigned> vertex_map;
    EEMP::SchreierTree spanning_tree;
    PermGroup schreier_generators;
  };

public:
  PartialPermInverseSemigroup();
  PartialPermInverseSemigroup(std::vector<PartialPerm> const &generators);

  void adjoin(std::vector<PartialPerm> const &generators, bool minimize = true);

  std::vector<PartialPerm> generators() const { return _generators; }
  bool empty() const { return _empty; }

  bool is_element(PartialPerm const &pperm) const;

private:
  void update_action_component(std::vector<PartialPerm> const &generators);
  void update_scc_representatives();

  bool _empty;

  std::vector<unsigned> _dom;
  std::vector<PartialPerm> _generators;

  std::vector<std::vector<unsigned>> _ac_im;
  std::unordered_map<std::vector<unsigned>, unsigned, VectorHash> _ac_im_ht;
  EEMP::SchreierTree _st_im;
  EEMP::OrbitGraph _og_im;

  std::vector<unsigned> _scc;
  std::vector<SccRepr> _scc_repr;

  std::vector<PartialPerm> _r_class_repr;
};

} // namespace cgtl

#endif // _GUARD_PARTIAL_PERM_INVERSE_SEMIGROUP_H
