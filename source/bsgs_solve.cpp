#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "bsgs.hpp"
#include "dbg.hpp"
#include "perm.hpp"
#include "perm_set.hpp"
#include "schreier_structure.hpp"

/**
 * @file bsgs_solve.cc
 * @brief Implements the solvable BSGS algorithm.
 *
 * @author Timo Nicolai
 */

namespace mpsym
{

namespace internal
{

void BSGS::solve(PermSet const &generators)
{
  DBG(DEBUG) << "Attempting to solve BSGS";

  unsigned iterations =
    static_cast<unsigned>(5.0 / 2.0 * std::log(degree()) / std::log(3.0));

  DBG(TRACE) << "Maximum number of iterations: " << iterations;

  for (auto const &gen : generators) {
    DBG(TRACE) << "Considering generator: " << gen;

    while (!strips_completely(gen)) {
      DBG(TRACE) << "Not in current BSGS";

      Perm w(gen);

      bool success = false;
      for (unsigned i = 0u; i < iterations; ++i) {
        DBG(TRACE) << "Iteration " << i;

        std::pair<Perm, Perm> conjugates;
        success = solve_s_normal_closure(generators, w, conjugates);
        if (success)
          break;

        Perm const &u(conjugates.first);
        Perm const &v(conjugates.second);
        DBG(TRACE) << "Conjugates are: " << u << " and " << v;

        w = ~u * ~v * u * v;
      }

      if (!success) {
        DBG(DEBUG) << "=> Failure";
        throw SolveError();
      }
    }
  }

  DBG(DEBUG) << "=> Success";
}

bool BSGS::solve_s_normal_closure(PermSet const &generators,
                                  Perm const &w,
                                  std::pair<Perm, Perm> &conjugates)
{
  DBG(TRACE) << "Begin calculating S-Normal Closure";

  BSGS const original_bsgs(*this);

  PermSet queue1 {w};
  PermSet queue2;

  for (auto i = 0u; i < queue1.size(); ++i) {
    Perm const g(queue1[i]);
    DBG(TRACE) << "Considering queue element: " << g;

    if (!strips_completely(g)) {
      DBG(TRACE) << "Not in current BSGS";

      for (auto const &h : queue2) {
        Perm tmp(~g * ~h * g * h);
        if (!original_bsgs.strips_completely(tmp)) {
          DBG(TRACE) << ~g << " * " << ~h << " * " << g << " * " << h
                     << " = " << tmp << " not in original BSGS";

          DBG(TRACE) << "=> Failure";
          DBG(TRACE) << "Finished calculating S-Normal Closure";

          conjugates.first = g;
          conjugates.second = h;

          return false;
        }
#ifndef NDEBUG
        else {
          DBG(TRACE) << ~g << " * " << ~h << " * " << g << " * " << h
                          << " = " << tmp << " in original BSGS";
        }
#endif
      }

      solve_adjoin_normalizing_generator(g);

      queue2.insert(g);

      DBG(TRACE) << "Updating queue:";
      for (auto const &gen : generators) {
        DBG(TRACE)
          << "  Appending: " << ~gen << " * " << g << " * " << gen
          << " = " << ~gen * g * gen;

        queue1.insert(~gen * g * gen);
      }
    }
#ifndef NDEBUG
    else
      DBG(TRACE) << "Already in current BSGS";
#endif
  }

  DBG(TRACE) << "=> Success";
  DBG(TRACE) << "Finished Calculating S-Normal Closure";
  return true;
}

void BSGS::solve_adjoin_normalizing_generator(Perm const &gen)
{
  DBG(TRACE) << "Begin adjoining normalizing generator";
  DBG(TRACE) << "Generator is: " << gen;

  unsigned i = 0u;
  Perm h(gen);

  while (!h.id()) {
    ++i;
    DBG(TRACE) << "Iteration " << i;

    if (i > base_size()) {
      for (unsigned j = 1u; j <= degree(); ++j) {
        if (h[j] != j) {
          extend_base(i);
          break;
        }
      }

      DBG(TRACE) << ">>> Updated base: " << _base << " <<<";
    }

    unsigned base_elem = base_point(i - 1u);

    throw std::logic_error("TODO: schreier structure initialization");

    DBG(TRACE)
      << "Considering h = " << h << " and b_" << i << " = " << base_elem
      << " (with orbit " << schreier_structure(i - 1u)->nodes() << ")";

    unsigned m = 1u;
    Perm h_m(h);
    unsigned tmp = h_m[base_elem];

    DBG(TRACE) << "h^1 = " << h_m;

    while (!schreier_structure(i - 1u)->contains(tmp)) {
      ++m;
      h_m *= h;
      tmp = h_m[base_elem];

      DBG(TRACE) << "h^" << m << " = " << h_m;
    }

    Perm u(schreier_structure(i - 1u)->transversal(tmp));
    DBG(TRACE) << "u = " << u;

    if (m > 1u) {
      DBG(TRACE) << "Enlarging:";

      for (unsigned j = 0u; j < i; ++j) { // TODO: avoid complete recomputation
        PermSet s_j(schreier_structure(j)->labels());
        s_j.insert(h);

        update_schreier_structure(j, s_j);

        DBG(TRACE) << "  S(" << j + 1u << ")" << " = " << s_j;
        DBG(TRACE) << "  O(" << j + 1u << ")" << " = " << schreier_structure(j)->nodes();
      }

      _strong_generators.insert(h);
      DBG(TRACE) << "  >>> Updated SGS: " << _strong_generators << " <<<";
    }

    h = h_m * ~u;
  }

  DBG(TRACE) << "Finished adjoining normalizing generator";
}

} // namespace internal

} // namespace mpsym
