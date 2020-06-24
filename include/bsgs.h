#ifndef _GUARD_BSGS_H
#define _GUARD_BSGS_H

#include <cassert>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include "orbits.h"
#include "perm.h"
#include "perm_set.h"
#include "schreier_generator_queue.h"
#include "schreier_structure.h"

namespace mpsym
{

class BSGSTransversalsBase
{
public:
  std::shared_ptr<SchreierStructure> schreier_structure(unsigned i) const
  { return _schreier_structures[i]; }

  void reserve_schreier_structure(
    unsigned i, unsigned root, unsigned degree);

  void update_schreier_structure(
    unsigned i, unsigned root, unsigned degree, PermSet const &generators);

  void insert_schreier_structure(
    unsigned i, unsigned root, unsigned degree, PermSet const &generators);

  void clear()
  { _schreier_structures.clear(); }

  virtual std::shared_ptr<SchreierStructure> make_schreier_structure(
    unsigned root, unsigned degree, PermSet const &generators) = 0;

private:
  std::vector<std::shared_ptr<SchreierStructure>> _schreier_structures;
};

template<typename T>
class BSGSTransversals : public BSGSTransversalsBase
{
  std::shared_ptr<SchreierStructure> make_schreier_structure(
    unsigned root, unsigned degree, PermSet const &generators) override
  { return std::make_shared<T>(degree, root, generators); }
};

struct BSGSOptions;

class BSGS
{
  friend std::ostream &operator<<(std::ostream &os, BSGS const &bsgs);

public:
  using order_type = boost::multiprecision::cpp_int;

  struct SolveError : public std::runtime_error
  {
    SolveError()
    : std::runtime_error("failed to solve BSGS")
    {}
  };

  explicit BSGS(unsigned degree = 1);

  BSGS(unsigned degree,
       PermSet const &generators,
       BSGSOptions const *options = nullptr);

  unsigned degree() const { return _degree; }
  order_type order() const;

  bool is_symmetric() const { return _is_symmetric; }
  bool is_alternating() const { return _is_alternating; }

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
  void construct_symmetric();
  void construct_alternating();
  void construct_unknown(PermSet const &generators, BSGSOptions const *options);

  // schreier sims initialization
  void schreier_sims(PermSet const &generators);

  void schreier_sims(std::vector<PermSet> &strong_generators,
                     std::vector<Orbit> &fundamental_orbits);

  void schreier_sims_random(PermSet const &generators,
                            BSGSOptions const *options);

  void schreier_sims_random(std::vector<PermSet> &strong_generators,
                            std::vector<Orbit> &fundamental_orbits,
                            BSGSOptions const *options);

  void schreier_sims_init(PermSet const &generators,
                          std::vector<PermSet> &strong_generators,
                          std::vector<Orbit> &fundamental_orbits);

  void schreier_sims_update_strong_gens(
    unsigned i,
    PermSet new_strong_generators,
    std::vector<PermSet> &strong_generators,
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

  std::shared_ptr<SchreierStructure> schreier_structure(unsigned i) const
  { return _transversals->schreier_structure(i); }

  void reserve_schreier_structure(unsigned i)
  {
    _transversals->reserve_schreier_structure(
      i, base_point(i), _degree);
  }

  void update_schreier_structure(unsigned i, PermSet const &generators)
  {
    _transversals->update_schreier_structure(
      i, base_point(i), _degree, generators);
  }

  void insert_schreier_structure(unsigned i, PermSet const &generators)
  {
    _transversals->insert_schreier_structure(
      i, base_point(i), _degree, generators);
  }

  unsigned _degree;

  std::vector<unsigned> _base;
  std::shared_ptr<BSGSTransversalsBase> _transversals;
  PermSet _strong_generators;

  bool _is_symmetric = false;
  bool _is_alternating = false;
};

std::ostream &operator<<(std::ostream &os, BSGS const &bsgs);

struct BSGSOptions
{
  enum class Construction {
    AUTO,
    SCHREIER_SIMS,
    SCHREIER_SIMS_RANDOM,
    SOLVE
  };

  enum class Transversals {
    EXPLICIT,
    SCHREIER_TREES,
    SHALLOW_SCHREIER_TREES
  };

  static BSGSOptions fill_defaults(BSGSOptions const *options)
  {
    static BSGSOptions default_options;
    return options ? *options : default_options;
  }

  Construction construction = Construction::AUTO;
  Transversals transversals = Transversals::EXPLICIT;

  bool check_altsym = true;
  bool reduce_gens = true;

  bool schreier_sims_random_guarantee = true;
  bool schreier_sims_random_use_known_order = true;
  BSGS::order_type schreier_sims_random_known_order = 0;
  int schreier_sims_random_retries = -1;
  unsigned schreier_sims_random_w = 100u;
};

} // namespace mpsym

#endif // _GUARD_BSGS_H
