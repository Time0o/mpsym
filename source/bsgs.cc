#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "schreier_sims.h"

namespace cgtl
{

namespace
{

void normalizing_generator(Perm const &gen, unsigned n, BSGS &bsgs)
{
  Dbg(Dbg::TRACE) << "=== BEGIN Adjoining normalizing generator";
  Dbg(Dbg::TRACE) << "Generator is: " << gen;

  unsigned i = 0u;
  Perm h(gen);

  while (!h.id()) {
    ++i;
    Dbg(Dbg::TRACE) << "== Iteration " << i;

    if (i > bsgs.base.size()) {
      for (unsigned j = 1u; j <= n; ++j) {
        if (h[j] != j) {
          bsgs.base.push_back(j);
          schreier_sims::SchreierTree st(n);
          st.create_root(j);
          bsgs.schreier_trees.push_back(st);
          break;
        }
      }

      Dbg(Dbg::TRACE) << ">>> Updated base: " << bsgs.base << " <<<";
    }

    unsigned base_elem = bsgs.base[i - 1u];

    schreier_sims::SchreierTree *schreier_tree = &bsgs.schreier_trees[i - 1u];

    Dbg(Dbg::TRACE)
      << "Considering h = " << h << " and b_" << i << " = " << base_elem
      << " (with orbit " << schreier_tree->nodes() << ")";

    unsigned m = 1u;
    Perm h_m(h);
    unsigned tmp = h_m[base_elem];

    Dbg(Dbg::TRACE) << "h^1 = " << h_m;

    while (!schreier_tree->contains(tmp)) {
      ++m;
      h_m *= h;
      tmp = h_m[base_elem];

      Dbg(Dbg::TRACE) << "h^" << m << " = " << h_m;
    }

    Perm u(schreier_tree->transversal(tmp));
    Dbg(Dbg::TRACE) << "u = " << u;

    if (m > 1u) {
      Dbg(Dbg::TRACE) << "Enlarging:";

      for (unsigned j = 0u; j < i; ++j) { // TODO: avoid complete recomputation
        std::vector<Perm> s_j(bsgs.schreier_trees[j].labels());
        s_j.push_back(h);

        schreier_sims::orbit(bsgs.base[j], s_j, &bsgs.schreier_trees[j]);

        Dbg(Dbg::TRACE) << "  S(" << j + 1u << ")" << " = " << s_j;
        Dbg(Dbg::TRACE) << "  O(" << j + 1u << ")" << " = "
                        << bsgs.schreier_trees[j].nodes();
      }

      bsgs.strong_generators.push_back(h);
      Dbg(Dbg::TRACE) << "  >>> Updated SGS: " << bsgs.strong_generators << " <<<";
    }

    h = h_m * ~u;
  }

  Dbg(Dbg::TRACE) << "=== END Adjoining normalizing generator";
}

bool s_normal_closure(
  Perm const &w, std::vector<Perm> const &generators, unsigned n,
  BSGS &bsgs, std::pair<Perm, Perm> &conjugates)
{
  Dbg(Dbg::TRACE) << "==== BEGIN Calculating S-Normal Closure";

  BSGS const original_bsgs(bsgs);

  std::vector<Perm> queue1 {w};
  std::vector<Perm> queue2;

  for (auto i = 0u; i < queue1.size(); ++i) {
    Perm const g(queue1[i]);
    Dbg(Dbg::TRACE) << "Considering queue element: " << g;

    if (!bsgs.contains(g)) {
      Dbg(Dbg::TRACE) << "=> Not in current BSGS";

      for (auto const &h : queue2) {
        Perm tmp(~g * ~h * g * h);
        if (!original_bsgs.contains(tmp)) {
          Dbg(Dbg::TRACE) << ~g << " * " << ~h << " * " << g << " * " << h
                          << " = " << tmp << " not in original BSGS";

          Dbg(Dbg::TRACE) << "==> Failure";
          Dbg(Dbg::TRACE) << "==== END Calculating S-Normal Closure";

          conjugates.first = g;
          conjugates.second = h;

          return false;
        }
#ifndef NDEBUG
        else {
          Dbg(Dbg::TRACE) << ~g << " * " << ~h << " * " << g << " * " << h
                          << " = " << tmp << " in original BSGS";
        }
#endif
      }

      normalizing_generator(g, n, bsgs);

      queue2.push_back(g);

      Dbg(Dbg::TRACE) << "=> Updating queue:";
      for (auto const &gen : generators) {
        Dbg(Dbg::TRACE)
          << "  Appending: " << ~gen << " * " << g << " * " << gen
          << " = " << ~gen * g * gen;

        queue1.push_back(~gen * g * gen);
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

}

bool BSGS::solve(std::vector<Perm> const &generators, BSGS &bsgs)
{
  assert(generators.size() > 0u);

  Dbg(Dbg::DBG) << "Attempting to solve BSGS for generators: " << generators;

  unsigned n = generators[0].degree();

  unsigned iterations =
    static_cast<unsigned>(5.0 / 2.0 * std::log(n) / std::log(3.0));

  Dbg(Dbg::TRACE) << "Maximum number of iterations: " << iterations;

  for (auto const &gen : generators) {
    Dbg(Dbg::TRACE) << "====== Considering generator: " << gen;

    while (!bsgs.contains(gen)) {
      Dbg(Dbg::TRACE) << "=> Not in current BSGS";

      Perm w(gen);

      bool success = false;
      for (unsigned i = 0u; i < iterations; ++i) {
        Dbg(Dbg::TRACE) << "===== Iteration " << i;

        std::pair<Perm, Perm> conjugates;
        success = s_normal_closure(w, generators, n, bsgs, conjugates);
        if (success)
          break;

        Perm const &u(conjugates.first);
        Perm const &v(conjugates.second);
        Dbg(Dbg::TRACE) << "=> Conjugates are: " << u << " and " << v;

        w = ~u * ~v * u * v;
      }

      if (!success) {
        Dbg(Dbg::DBG) << "==> Failure";
        return false;
      }
    }
  }

  Dbg(Dbg::DBG) << "==> Success";
  return true;
}

} // namespace cgtl
