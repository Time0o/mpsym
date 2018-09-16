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

std::vector<Perm> normalizing_generator(Perm const &gen, BSGS &bsgs)
{
  unsigned n = bsgs.strong_generators[0].degree();

  // returned list of newly appended strong generators
  std::vector<Perm> res;

  // main loop
  unsigned i = 0u;
  Perm h(gen);

  while (!h.id()) {
    if (i >= bsgs.base.size()) {
      for (unsigned j = 1u; j <= n; ++j) {
        if (h[j] != j) {
          bsgs.base.push_back(j);
          bsgs.schreier_trees.push_back(schreier_sims::SchreierTree(n));
        }
      }

      Perm h_m(h);
      unsigned tmp = h_m[bsgs.base[i]];

      bool enlarge = false;
      while (!bsgs.schreier_trees[i].contains(tmp)) {
        h_m *= h_m;
        tmp = h_m[bsgs.base[i]];

        if (!enlarge)
          enlarge = true;
      }

      Perm u(bsgs.schreier_trees[i].transversal(tmp));

      if (enlarge) {
        std::vector<Perm> s_i(bsgs.schreier_trees[i].labels());
        s_i.push_back(h);

        schreier_sims::orbit(bsgs.base[i], s_i,  &bsgs.schreier_trees[i]);
        res.push_back(h);
      }

      h = h_m * ~u;
      ++i;
    }
  }

  return res;
}

SNormalClosure s_normal_closure(
  Perm const &w, std::vector<Perm> const &generators, BSGS &bsgs)
{
  std::vector<Perm> res;

  BSGS original_bsgs(bsgs);

  std::vector<Perm> queue1 {w};
  std::vector<Perm> queue2;

  for (auto i = 0u; i < queue1.size(); ++i) {
    Perm const &g(queue1[i]);

    if (!bsgs.contains(g)) {
      for (auto const &h : queue2) {
        if (!original_bsgs.contains(~g * ~h * g * h))
          return SNormalClosure(g, h);
      }

      auto tmp(normalizing_generator(g, bsgs));
      res.insert(res.end(), tmp.begin(), tmp.end());

      queue2.push_back(g);

      for (auto const &gen : generators)
        queue1.push_back(~gen * g * gen);
    }
  }

  return SNormalClosure(res);
}

}

bool BSGS::solve(std::vector<unsigned> const &partial_base,
                 std::vector<Perm> const &generators, BSGS &bsgs)
{
  unsigned n = generators[0].degree();

  bsgs.base = partial_base;
  bsgs.strong_generators = std::vector<Perm>();

  bsgs.schreier_trees =
    std::vector<schreier_sims::SchreierTree>(bsgs.base.size(), n);

  for (auto i = 0u; i < bsgs.base.size(); ++i)
    bsgs.schreier_trees[i].create_root(bsgs.base[i]);

  unsigned iterations =
    static_cast<unsigned>(5.0 / 2.0 * std::log(n) / std::log(3));

  for (auto const &gen : generators) {
    while (!bsgs.contains(gen)) {
      Perm w(gen);

      bool success = false;
      for (unsigned i = 0u; i < iterations; ++i) {
        SNormalClosure sncl(s_normal_closure(w, generators, bsgs));
        if (sncl.success) {
          success = true;
          break;
        }

        Perm const &u(sncl.conjugates.first);
        Perm const &v(sncl.conjugates.second);

        w = ~u * ~v * u * v;
      }

      if (!success) {
        return false;
      }
    }
  }

  return true;
}

} // namespace cgtl
