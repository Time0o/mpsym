#include <cassert>
#include <vector>

#include "dbg.h"
#include "perm.h"
#include "schreier_sims.h"
#include "solvable_bsgs.h"

namespace cgtl
{

std::vector<Perm> normalizing_generator(
  Perm const &gen, std::vector<unsigned> &base, std::vector<Perm> &generators,
  std::vector<schreier_sims::SchreierTree> &sts)
{
#ifndef NDEBUG
  for (auto const &_gen : generators)
    assert(_gen.degree() == generators[0].degree());
#endif

  unsigned n = generators[0].degree();

  // TODO: assert that gen normalizes?

  // returned list of newly appended strong generators
  std::vector<Perm> res;

  // main loop
  unsigned i = 0u;
  Perm h(gen);

  while (!h.id()) {
    if (i >= base.size()) {
      for (unsigned j = 1u; j <= n; ++j) {
        if (h[j] != j) {
          base.push_back(j);
          sts.push_back(schreier_sims::SchreierTree(n));
        }
      }

      Perm h_m(h);
      unsigned tmp = h_m[base[i]];

      bool enlarge = false;
      while (!sts[i].contains(tmp)) {
        h_m *= h_m;
        tmp = h_m[base[i]];

        if (!enlarge)
          enlarge = true;
      }

      Perm u(sts[i].transversal(tmp));

      if (enlarge) {
        std::vector<Perm> s_i(sts[i].labels());
        s_i.push_back(h);

        schreier_sims::orbit(base[i], s_i,  &sts[i]);
        res.push_back(h);
      }

      h = h_m * ~u;
      ++i;
    }
  }

  return res;
}

} // namespace cgtl
