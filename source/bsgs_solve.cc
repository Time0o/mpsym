#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "perm_set.h"
#include "schreier_structure.h"

/**
 * @file bsgs_solve.cc
 * @brief Implements the solvable BSGS algorithm.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

void BSGS::solve(PermSet const &generators)
{
  Dbg(Dbg::DBG) << "Attempting to solve BSGS for generators: " << generators;

  unsigned iterations =
    static_cast<unsigned>(5.0 / 2.0 * std::log(degree()) / std::log(3.0));

  Dbg(Dbg::TRACE) << "Maximum number of iterations: " << iterations;

  for (auto const &gen : generators) {
    Dbg(Dbg::TRACE) << "====== Considering generator: " << gen;

    while (!strips_completely(gen)) {
      Dbg(Dbg::TRACE) << "=> Not in current BSGS";

      Perm w(gen);

      bool success = false;
      for (unsigned i = 0u; i < iterations; ++i) {
        Dbg(Dbg::TRACE) << "===== Iteration " << i;

        std::pair<Perm, Perm> conjugates;
        success = solve_s_normal_closure(generators, w, &conjugates);
        if (success)
          break;

        Perm const &u(conjugates.first);
        Perm const &v(conjugates.second);
        Dbg(Dbg::TRACE) << "=> Conjugates are: " << u << " and " << v;

        w = ~u * ~v * u * v;
      }

      if (!success) {
        Dbg(Dbg::DBG) << "==> Failure";
        throw SolveError();
      }
    }
  }

  Dbg(Dbg::DBG) << "==> Success";
}

bool BSGS::solve_s_normal_closure(PermSet const &generators,
                                  Perm const &w,
                                  std::pair<Perm, Perm> *conjugates)
{
  Dbg(Dbg::TRACE) << "==== BEGIN Calculating S-Normal Closure";

  BSGS const original_bsgs(*this);

  PermSet queue1 {w};
  PermSet queue2;

  for (auto i = 0u; i < queue1.size(); ++i) {
    Perm const g(queue1[i]);
    Dbg(Dbg::TRACE) << "Considering queue element: " << g;

    if (!strips_completely(g)) {
      Dbg(Dbg::TRACE) << "=> Not in current BSGS";

      for (auto const &h : queue2) {
        Perm tmp(~g * ~h * g * h);
        if (!original_bsgs.strips_completely(tmp)) {
          Dbg(Dbg::TRACE) << ~g << " * " << ~h << " * " << g << " * " << h
                          << " = " << tmp << " not in original BSGS";

          Dbg(Dbg::TRACE) << "==> Failure";
          Dbg(Dbg::TRACE) << "==== END Calculating S-Normal Closure";

          conjugates->first = g;
          conjugates->second = h;

          return false;
        }
#ifndef NDEBUG
        else {
          Dbg(Dbg::TRACE) << ~g << " * " << ~h << " * " << g << " * " << h
                          << " = " << tmp << " in original BSGS";
        }
#endif
      }

      solve_adjoin_normalizing_generator(g);

      queue2.add(g);

      Dbg(Dbg::TRACE) << "=> Updating queue:";
      for (auto const &gen : generators) {
        Dbg(Dbg::TRACE)
          << "  Appending: " << ~gen << " * " << g << " * " << gen
          << " = " << ~gen * g * gen;

        queue1.add(~gen * g * gen);
      }
    }
#ifndef NDEBUG
    else
      Dbg(Dbg::TRACE) << "=> Already in current BSGS";
#endif
  }

  Dbg(Dbg::TRACE) << "==> Success";
  Dbg(Dbg::TRACE) << "==== END Calculating S-Normal Closure";
  return true;
}

void BSGS::solve_adjoin_normalizing_generator(Perm const &gen)
{
  Dbg(Dbg::TRACE) << "=== BEGIN Adjoining normalizing generator";
  Dbg(Dbg::TRACE) << "Generator is: " << gen;

  unsigned i = 0u;
  Perm h(gen);

  while (!h.id()) {
    ++i;
    Dbg(Dbg::TRACE) << "== Iteration " << i;

    if (i > base.size()) {
      for (unsigned j = 1u; j <= degree(); ++j) {
        if (h[j] != j) {
          extend_base(i);
          break;
        }
      }

      Dbg(Dbg::TRACE) << ">>> Updated base: " << base << " <<<";
    }

    unsigned base_elem = base[i - 1u];

    auto schreier_structure(schreier_structures[i - 1u]);

    Dbg(Dbg::TRACE)
      << "Considering h = " << h << " and b_" << i << " = " << base_elem
      << " (with orbit " << schreier_structure->nodes() << ")";

    unsigned m = 1u;
    Perm h_m(h);
    unsigned tmp = h_m[base_elem];

    Dbg(Dbg::TRACE) << "h^1 = " << h_m;

    while (!schreier_structure->contains(tmp)) {
      ++m;
      h_m *= h;
      tmp = h_m[base_elem];

      Dbg(Dbg::TRACE) << "h^" << m << " = " << h_m;
    }

    Perm u(schreier_structure->transversal(tmp));
    Dbg(Dbg::TRACE) << "u = " << u;

    if (m > 1u) {
      Dbg(Dbg::TRACE) << "Enlarging:";

      for (unsigned j = 0u; j < i; ++j) { // TODO: avoid complete recomputation
        PermSet s_j(schreier_structures[j]->labels());
        s_j.add(h);

        update_schreier_structure(j, s_j.vect());

        Dbg(Dbg::TRACE) << "  S(" << j + 1u << ")" << " = " << s_j;
        Dbg(Dbg::TRACE) << "  O(" << j + 1u << ")" << " = "
                        << schreier_structures[j]->nodes();
      }

      strong_generators.push_back(h);
      Dbg(Dbg::TRACE) << "  >>> Updated SGS: " << strong_generators << " <<<";
    }

    h = h_m * ~u;
  }

  Dbg(Dbg::TRACE) << "=== END Adjoining normalizing generator";
}

} // namespace cgtl
