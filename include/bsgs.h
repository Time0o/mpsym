#ifndef _GUARD_BSGS_H
#define _GUARD_BSGS_H

#include <cassert>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "orbits.h"
#include "perm.h"
#include "perm_set.h"
#include "schreier_generator_queue.h"
#include "schreier_structure.h"

namespace cgtl
{

class BSGS
{
  friend std::ostream &operator<<(std::ostream &os, BSGS const &bsgs);

public:
  struct SolveError : public std::runtime_error
  {
    SolveError()
    : std::runtime_error("failed to solve BSGS")
    {}
  };

  enum class Construction {
    SCHREIER_SIMS,
    SCHREIER_SIMS_RANDOM,
    SOLVE,
    AUTO
  };

  enum class Transversals {
    EXPLICIT,
    SCHREIER_TREES,
    SHALLOW_SCHREIER_TREES,
    AUTO
  };

  BSGS(unsigned degree,
       PermSet const &generators,
       Construction construction = Construction::AUTO,
       Transversals transversals = Transversals::AUTO);

  unsigned degree() const { return _degree; }

  bool base_empty() const { return _base.empty(); }
  unsigned base_size() const { return _base.size(); }
  unsigned base_point(unsigned i) const { return _base[i]; }
  void base_change(std::vector<unsigned> prefix);

  PermSet strong_generators() const { return _strong_generators; }
  PermSet strong_generators(unsigned i) const;

  Orbit orbit(unsigned i) const;
  Perm transversal(unsigned i, unsigned o) const;
  PermSet transversals(unsigned i) const;
  PermSet stabilizers(unsigned i) const;

  std::pair<Perm, unsigned> strip(Perm const &perm, unsigned offs = 0) const;
  bool strips_completely(Perm const &perm) const;

private:
  // schreier sims initialization
  void schreier_sims(PermSet const &generators);
  void schreier_sims_random(PermSet const &generators, unsigned w = 10);

  void schreier_sims_init(std::vector<PermSet> &strong_generators,
                          std::vector<Orbit> &fundamental_orbits);

  void schreier_sims_finish();

  // solvable BSGS initialization
  void solve(PermSet const &generators);

  bool solve_s_normal_closure(PermSet const &generators,
                              Perm const &w,
                              std::pair<Perm, Perm> &conjugates);

  void solve_adjoin_normalizing_generator(Perm const &gen);

  // generator reduction
  void reduce_gens();

  std::unordered_set<Perm> reduce_gens_set_difference(
    std::unordered_set<Perm> const &lhs,
    std::unordered_set<Perm> const &rhs,
    std::unordered_set<Perm> const &base) const;

  // base change
  void swap_base_points(unsigned i);
  void transpose_base_point(unsigned i, unsigned j);
  unsigned insert_redundant_base_point(unsigned bp, unsigned i_min);
  void conjugate(Perm const &conj);

  // convenience methods
  void extend_base(unsigned bp);
  void extend_base(unsigned bp, unsigned i);

  void update_schreier_structure(unsigned i, PermSet const &strong_generators);
  void insert_schreier_structure(unsigned i, PermSet const &strong_generators);

  unsigned _degree;
  Transversals _transversals;
  std::vector<unsigned> _base;
  std::vector<std::shared_ptr<SchreierStructure>> _schreier_structures;
  PermSet _strong_generators;
};

std::ostream &operator<<(std::ostream &os, BSGS const &bsgs);

} // namespace cgtl

#endif // _GUARD_BSGS_H
