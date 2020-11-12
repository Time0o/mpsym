#ifndef GUARD_PERM_GROUP_H
#define GUARD_PERM_GROUP_H

#include <cassert>
#include <map>
#include <tuple>
#include <type_traits>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include "bsgs.hpp"
#include "iterator.hpp"
#include "perm.hpp"
#include "perm_set.hpp"
#include "timeout.hpp"

namespace mpsym
{

namespace internal
{

class Orbit;
class OrbitPartition;

class BlockSystem;

class PermGroup
{
  friend std::ostream &operator<<(std::ostream &os, PermGroup const &pg);

public:
  class const_iterator : public ForwardIterator<const_iterator, Perm const>
  {
  public:
    const_iterator() : _end(true) {};
    const_iterator(PermGroup const &pg);

    bool operator==(const_iterator const &rhs) const override;

    PermSet const &factors() const
    { return _current_factors; }

  private:
    reference current() override;
    void next() override;

    std::vector<unsigned> _state;
    bool _trivial;
    bool _end;

    std::vector<PermSet> _transversals;
    Perm _current;
    bool _current_valid;
    PermSet _current_factors;
  };

  explicit PermGroup(unsigned degree = 1)
  : _bsgs(degree),
    _order(1)
  {}

  explicit PermGroup(BSGS const &bsgs)
  : _bsgs(bsgs),
    _order(bsgs.order())
  {}

  PermGroup(unsigned degree, PermSet const &generators);

  static PermGroup symmetric(unsigned degree);
  static PermGroup cyclic(unsigned degree);
  static PermGroup dihedral(unsigned degree);

  template<typename IT>
  static PermGroup direct_product(IT first,
                                  IT last,
                                  BSGSOptions const *bsgs_options_ = nullptr,
                                  timeout::flag aborted = timeout::unset())
  {
    assert(std::distance(first, last) > 0);

    // direct product degree
    unsigned dp_degree = 0u;

    for (auto it = first; it != last; ++it)
      dp_degree += it->degree();

    // direct product order
    auto dp_order(direct_product_order(first, last));

    // direct product generators
    unsigned d = 0u;

    PermSet dp_generators;
    for (auto it = first; it != last; ++it) {
      for (Perm const &perm : it->generators())
        dp_generators.insert(perm.shifted(d).extended(dp_degree));

      d += it->degree();
    }

    // construct direct product
    auto bsgs_options(BSGSOptions::fill_defaults(bsgs_options_));
    bsgs_options.schreier_sims_random_known_order = dp_order;

    return PermGroup(BSGS(dp_degree, dp_generators, &bsgs_options, aborted));
  }

  template<typename IT>
  static BSGS::order_type direct_product_order(IT first, IT last)
  {
    BSGS::order_type dp_order = 1u;

    for (auto it = first; it != last; ++it)
      dp_order *= it->order();

    return dp_order;
  }

  static PermGroup wreath_product(PermGroup const &lhs,
                                  PermGroup const &rhs,
                                  BSGSOptions const *bsgs_options = nullptr,
                                  timeout::flag aborted = timeout::unset());

  static BSGS::order_type wreath_product_order(PermGroup const &lhs,
                                               PermGroup const &rhs);

  bool operator==(PermGroup const &rhs) const;
  bool operator!=(PermGroup const &rhs) const;

  const_iterator begin() const { return const_iterator(*this); }
  const_iterator end() const { return const_iterator(); }

  PermSet generators() const { return _bsgs.strong_generators(); }

  BSGS &bsgs() { return _bsgs; }
  BSGS const &bsgs() const { return _bsgs; }

  unsigned degree() const { return _bsgs.degree(); }
  BSGS::order_type order() const { return _order; }

  unsigned smallest_moved_point() const
  { return generators().smallest_moved_point(); }

  unsigned largest_moved_point() const
  { return generators().largest_moved_point(); }

  bool is_trivial() const { return _bsgs.base_empty(); }
  bool is_symmetric() const;
  bool is_shifted_symmetric() const;
  bool is_transitive() const;

  bool contains_element(Perm const &perm) const;
  Perm random_element() const;

  std::vector<PermGroup> disjoint_decomposition(
    bool complete = true, bool disjoint_orbit_optimization = false) const;

  std::vector<PermGroup> wreath_decomposition() const;

private:
  static boost::multiprecision::cpp_int symmetric_order(unsigned deg)
  {
    boost::multiprecision::cpp_int ret(1);
    for (unsigned i = deg; i > 0u; --i)
      ret *= i;

    return ret;
  }

  // complete disjoint decomposition
  bool disjoint_decomp_orbits_dependent(
    Orbit const &orbit1,
    Orbit const &orbit2) const;

  void disjoint_decomp_generate_dependency_classes(
    OrbitPartition &orbits) const;

  static bool disjoint_decomp_restricted_subgroups(
    OrbitPartition const &orbit_split,
    PermGroup const &perm_group,
    std::pair<PermGroup, PermGroup> &restricted_subgroups);

  static std::vector<PermGroup> disjoint_decomp_join_results(
    std::vector<PermGroup> const &res1,
    std::vector<PermGroup> const &res2);

  static std::vector<PermGroup> disjoint_decomp_complete_recursive(
    OrbitPartition const &orbits,
    PermGroup const &perm_group);

  std::vector<PermGroup> disjoint_decomp_complete(
    bool disjoint_orbit_optimization = true) const;

  // incomplete disjoint decomposition
  struct MovedSet : public std::vector<unsigned>
  {
    void init(Perm const &perm);
    bool equivalent(MovedSet const &other) const;
    void extend(MovedSet const &other);
  };

  struct EquivalenceClass
  {
    EquivalenceClass(Perm const &init, MovedSet const &moved)
    : generators({init}),
      moved(moved),
      merged(false)
    {}

    PermSet generators;
    MovedSet moved;
    bool merged;
  };

  std::vector<EquivalenceClass> disjoint_decomp_find_equivalence_classes() const;

  void disjoint_decomp_merge_equivalence_classes(
    std::vector<EquivalenceClass> &equivalence_classes) const;

  std::vector<PermGroup> disjoint_decomp_incomplete() const;

  // wreath decomposition
  std::vector<PermGroup> wreath_decomp_find_stabilizers(
    BlockSystem const &block_system,
    PermGroup const &block_permuter) const;

  PermSet wreath_decomp_construct_block_permuter_image(
    BlockSystem const &block_system,
    PermGroup const &block_permuter) const;

  bool wreath_decomp_reconstruct_block_permuter(
    BlockSystem const &block_system,
    PermGroup const &block_permuter,
    PermSet const &block_permuter_image) const;

  BSGS _bsgs;
  BSGS::order_type _order;
};

std::ostream &operator<<(std::ostream &os, PermGroup const &pg);

} // namespace internal

} // namespace mpsym

#endif // GUARD_PERM_GROUP_H
