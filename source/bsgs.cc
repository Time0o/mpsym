#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "schreier_sims.h"
#include "schreier_structure.h"

namespace cgtl
{

bool BSGS::contains(Perm const &perm) const {
  auto strip_result(schreier_sims::strip(perm, base, schreier_structures));
  return strip_result.first.id() && strip_result.second == base.size() + 1u;
}

std::vector<unsigned> BSGS::orbit(unsigned i) const {
  return schreier_structures[i]->nodes();
}

Perm BSGS::transversal(unsigned i, unsigned o) const {
  return schreier_structures[i]->transversal(o);
}

std::vector<Perm> BSGS::transversals(unsigned i) const {
  std::vector<Perm> transversals;
  for (unsigned o : orbit(i))
    transversals.push_back(schreier_structures[i]->transversal(o));

  return transversals;
}

std::vector<Perm> BSGS::stabilizers(unsigned i) const {
  return schreier_structures[i]->labels();
}

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
          SchreierTree st(n);
          st.create_root(j);
          bsgs.schreier_structures.push_back(std::make_shared<SchreierTree>(st));
          break;
        }
      }

      Dbg(Dbg::TRACE) << ">>> Updated base: " << bsgs.base << " <<<";
    }

    unsigned base_elem = bsgs.base[i - 1u];

    std::shared_ptr<SchreierStructure> schreier_structure(
      bsgs.schreier_structures[i - 1u]);

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
        std::vector<Perm> s_j(bsgs.schreier_structures[j]->labels());
        s_j.push_back(h);

        schreier_sims::orbit(bsgs.base[j], s_j, bsgs.schreier_structures[j]);

        Dbg(Dbg::TRACE) << "  S(" << j + 1u << ")" << " = " << s_j;
        Dbg(Dbg::TRACE) << "  O(" << j + 1u << ")" << " = "
                        << bsgs.schreier_structures[j]->nodes();
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

BSGS BSGS::solve(std::vector<Perm> const &generators)
{
  assert(generators.size() > 0u);

  Dbg(Dbg::DBG) << "Attempting to solve BSGS for generators: " << generators;

  BSGS bsgs;

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
        bsgs.base.clear();
        bsgs.strong_generators.clear();
        bsgs.schreier_structures.clear();
        return bsgs;
      }
    }
  }

  Dbg(Dbg::DBG) << "==> Success";
  return bsgs;
}

} // namespace cgtl
