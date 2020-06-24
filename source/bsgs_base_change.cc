#include <algorithm>
#include <cassert>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"

namespace mpsym
{

namespace internal
{

void BSGS::base_change(std::vector<unsigned> prefix)
{
  DBG(DEBUG) << "Appending prefix " << prefix << " to base " << _base;

  Perm conj(degree());
  Perm conj_inv(degree());

  for (auto i = 0u; i < prefix.size(); ++i) {
    unsigned target = conj_inv[prefix[i]];

    if (i >= base_size()) {
      DBG(TRACE) << "Prefix point: " << prefix[i];

      insert_redundant_base_point(target, i);
      DBG(TRACE) << "Appended " << target << " to base: " << _base;

      continue;
    }

    DBG(TRACE) << "Base/prefix points: "
               << base_point(i) << "/" << prefix[i] << "(" << target << ")";

    if (base_point(i) == target)
      continue;

    if (schreier_structure(i)->contains(target)) {
      // update the conjugation permutation so that it will correctly conjugate
      // all base points until position i to the corresponsing prefix base points
      Perm transv(schreier_structure(i)->transversal(target));

      DBG(TRACE) << target << " in O(" << i + 1u << ") = " << orbit(i)
                 << " (transversal is " << transv << ")";

      conj = transv * conj;
      conj_inv = ~conj;

      DBG(TRACE) << "Updated conjugating permutation: " << conj;

    } else {
      DBG(TRACE) << target << " not in O(" << i + 1u << ") = " << orbit(i);

      unsigned j = insert_redundant_base_point(target, i);
      DBG(TRACE) << "Inserted " << target << " into base: " << _base;

      transpose_base_point(j, i);
      DBG(TRACE) << "Base after transposition: " << _base;
    }
  }

  DBG(TRACE) << "Final conjugating permutation: " << conj;

  conjugate(conj);

  DBG(DEBUG) << "Base after conjugation:" << _base;

  assert(std::equal(prefix.begin(), prefix.end(), _base.begin()));
}

void BSGS::swap_base_points(unsigned i)
{
  DBG(TRACE) << "Swapping base points " << i + 1u << " and " << i + 2u;

  assert(i < base_size() - 1u);

  // swap base point values
  DBG(TRACE) << "Previous base: " << _base;
  std::swap(_base[i], _base[i + 1u]);
  DBG(TRACE) << "New base: " << _base;

  // recompute schreier structures for base points i and i+1
  auto sgi(stabilizers(i));
  auto oi(orbit(i));

  update_schreier_structure(i, sgi);

  auto sgi1(strong_generators(i + 1u));
  auto oi1(orbit(i + 1u));

  update_schreier_structure(i + 1u, sgi1);

  // final size of the fundamental orbit O_(i+1)
  auto oi1_desired_size = (oi.size() * oi1.size()) / orbit(i).size();

  DBG(TRACE) << "Desired size of O(" << i + 1u << ") is " << oi1_desired_size;

  // iterate over schreier generators
  SchreierGeneratorQueue schreier_generator_queue;

  sgi = stabilizers(i);
  oi = orbit(i);

  schreier_generator_queue.update(sgi, oi, schreier_structure(i));

  for (Perm const &perm : schreier_generator_queue) {
    DBG(TRACE) << "Schreier Generator: " << perm;

    if (!schreier_structure(i + 1)->contains(perm[base_point(i + 1u)])) {
      DBG(TRACE) << "Updating strong generators:";

      // extend strong generators
      sgi1.insert(perm);
      update_schreier_structure(i + 1u, sgi1);

      DBG(TRACE) << "S(" << i + 1u << ") = " << stabilizers(i + 1u);
      DBG(TRACE) << "O(" << i + 1u << ") = " << orbit(i + 1u);

      if (orbit(i + 1u).size() >= oi1_desired_size)
        break;
    }
  }

  assert(orbit(i + 1u).size() >= oi1_desired_size);

  DBG(TRACE) << "Final size of O(" << i + 1u << ") is " << orbit(i + 1u).size();

  // eliminate duplicate strong generators
  _strong_generators.insert(sgi1.begin(), sgi1.end());
  _strong_generators.make_unique();
}

void BSGS::transpose_base_point(unsigned i, unsigned j)
{
   // swap base points until the base point at position i has moved to position j
   while (i > j) {
     swap_base_points(i - 1u);
     --i;
   }
}

unsigned BSGS::insert_redundant_base_point(unsigned bp, unsigned i_min)
{
  unsigned i = std::min(i_min + 1u, static_cast<unsigned>(base_size()));

  while (i < base_size()) {
    if (base_point(i) == bp)
      return i;

    bool stabilized = true;

    for (Perm const &stab : stabilizers(i - 1u)) {
       if (stab[bp] != bp) {
          stabilized = false;
          break;
      }
    }

    if (stabilized)
      break;

    ++i;
  }

  bool reuse_stabilizers = i < base_size();

  // insert new base point
  extend_base(bp, i);

  // compute schreier structure for new base point
  insert_schreier_structure(
    i, reuse_stabilizers ? stabilizers(i - 1u) : strong_generators(i));

  return i;
}

void BSGS::conjugate(Perm const &conj)
{
  // conjugate base
  for (unsigned &b : _base)
    b = conj[b];

  // conjugate strong generating set
  for (Perm &sg : _strong_generators)
    sg = ~conj * sg * conj;

  // update schreier structures
  for (unsigned i = 0u; i < base_size(); ++i)
    update_schreier_structure(i, strong_generators(i));
}

} // namespace internal

} // namespace mpsym
