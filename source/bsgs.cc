#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "schreier_sims.h"
#include "solvable_bsgs.h"

namespace cgtl
{

namespace
{

struct SNormalClosure {
  SNormalClosure(std::vector<Perm> const &result)
    : success(true), result(result), conjugates() {}

  SNormalClosure(Perm const &u, Perm const &v)
    : success(false), result(), conjugates(std::make_pair(u, v)) {}

  bool success;
  std::vector<Perm> result;
  std::pair<Perm, Perm> conjugates;
};

std::vector<Perm> normalizing_generator(Perm const &gen, unsigned n, BSGS &bsgs)
{
  Dbg(Dbg::TRACE) << "= Adjoining normalizing generator: " << gen;

  std::vector<Perm> res;

  unsigned i = 0u;
  Perm h(gen);

  while (!h.id()) {
    if (++i >= bsgs.base.size()) {
      Dbg(Dbg::TRACE) << "i = " << i;

      for (unsigned j = 1u; j <= n; ++j) {
        if (h[j] != j) {
          bsgs.base.push_back(j);
          schreier_sims::SchreierTree st(n);
          st.create_root(j);
          bsgs.schreier_trees.push_back(st);
          break;
        }
      }

      Dbg(Dbg::TRACE) << "Updated base: " << bsgs.base;
    }

    unsigned base_elem = bsgs.base[i - 1u];

    Dbg(Dbg::TRACE) << "Considering permutation: " << h;
    Dbg(Dbg::TRACE) << "Considering base element: " << base_elem;

    schreier_sims::SchreierTree *schreier_tree = &bsgs.schreier_trees[i - 1u];

    Perm h_m(h);
    Dbg(Dbg::TRACE) << "h_m = " << h_m;

    unsigned tmp = h_m[base_elem];

    bool enlarge = false;
    while (!schreier_tree->contains(tmp)) {
      h_m = h_m * h_m;
      Dbg(Dbg::TRACE) << "h_m = " << h_m;

      tmp = h_m[base_elem];

      if (!enlarge)
        enlarge = true;
    }

    Perm u(schreier_tree->transversal(tmp));
    Dbg(Dbg::TRACE) << "u = " << h_m;

    if (enlarge) {
      Dbg(Dbg::TRACE) << "Enlarging...";

      std::vector<Perm> s_i(schreier_tree->labels());
      s_i.push_back(h);

      schreier_sims::orbit(base_elem, s_i, schreier_tree);
      res.push_back(h);
    }

    h = h_m * ~u;
  }

  Dbg(Dbg::TRACE) << "==> Result is: " << res;
  Dbg(Dbg::TRACE) << "=";
  return res;
}

SNormalClosure s_normal_closure(
  Perm const &w, std::vector<Perm> const &generators, unsigned n, BSGS &bsgs)
{
  Dbg(Dbg::TRACE) << "== Calculating S-Normal Closure";

  std::vector<Perm> res;

  BSGS original_bsgs(bsgs);

  std::vector<Perm> queue1 {w};
  std::vector<Perm> queue2;

  for (auto i = 0u; i < queue1.size(); ++i) {
    Perm const g(queue1[i]);
    Dbg(Dbg::TRACE) << "Considering queue element: " << g;

    if (!bsgs.contains(g)) {
      Dbg(Dbg::TRACE) << "=> Not in current BSGS";

      for (auto const &h : queue2) {
        Dbg(Dbg::TRACE) << "Checking: " << h;

        if (!original_bsgs.contains(~g * ~h * g * h)) {
          Dbg(Dbg::TRACE) << "==> Not in original BSGS, failure";
          Dbg(Dbg::TRACE) << "==";
          return SNormalClosure(g, h);
        }
      }

      auto tmp(normalizing_generator(g, n, bsgs));

      res.insert(res.end(), tmp.begin(), tmp.end());
      Dbg(Dbg::TRACE) << "=> Updated result is: " << res;

      queue2.push_back(g);

      for (auto const &gen : generators) {
        Dbg(Dbg::TRACE) << "Appending to queue: " << ~gen * g * gen;
        queue1.push_back(~gen * g * gen);
      }
    }
  }

  Dbg(Dbg::TRACE) << "==> Success, result is: " << res;
  Dbg(Dbg::TRACE) << "==";
  return SNormalClosure(res);
}

}

bool BSGS::solve(std::vector<Perm> const &generators, BSGS &bsgs)
{
  assert(generators.size() > 0u);

  Dbg(Dbg::DBG) << "Attempting to solve BSGS for generators: " << generators;

  unsigned n = generators[0].degree();

  unsigned iterations =
    static_cast<unsigned>(5.0 / 2.0 * std::log(n) / std::log(3));

  Dbg(Dbg::TRACE) << "Maximum number of iterations: " << iterations;

  for (auto const &gen : generators) {
    Dbg(Dbg::TRACE) << "==== Considering generator: " << gen;

    while (!bsgs.contains(gen)) {
      Dbg(Dbg::TRACE) << "=> Not in current BSGS";

      Perm w(gen);
      Dbg(Dbg::TRACE) << "=> w := " << w;

      bool success = false;
      for (unsigned i = 0u; i < iterations; ++i) {
        Dbg(Dbg::TRACE) << "=== Iteration " << i;

        SNormalClosure sncl(s_normal_closure(w, generators, n, bsgs));
        if (sncl.success) {
          success = true;
          break;
        }

        Perm const &u(sncl.conjugates.first);
        Perm const &v(sncl.conjugates.second);
        Dbg(Dbg::TRACE) << "=> Conjugates are: " << u << " and " << v;

        w = ~u * ~v * u * v;
        Dbg(Dbg::TRACE) << "=> w := " << w;
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
