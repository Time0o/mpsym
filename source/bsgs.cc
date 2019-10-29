#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <ostream>
#include <unordered_set>
#include <utility>
#include <vector>

#include "bsgs.h"
#include "dbg.h"
#include "perm.h"
#include "schreier_sims.h"
#include "schreier_structure.h"

namespace cgtl
{

std::vector<unsigned> BSGS::orbit(unsigned i) const
{
  return schreier_structures[i]->nodes();
}

Perm BSGS::transversal(unsigned i, unsigned o) const
{
  return schreier_structures[i]->transversal(o);
}

std::vector<Perm> BSGS::transversals(unsigned i) const
{
  std::vector<Perm> transversals;
  for (unsigned o : orbit(i))
    transversals.push_back(schreier_structures[i]->transversal(o));

  return transversals;
}

std::vector<Perm> BSGS::stabilizers(unsigned i) const
{
  return schreier_structures[i]->labels();
}

std::pair<Perm, unsigned> BSGS::strip(Perm const &perm, unsigned offs) const
{
  Perm result(perm);

  for (unsigned i = offs; i < base.size(); ++i) {
    unsigned beta = result[base[i]];
    if (!schreier_structures[i]->contains(beta))
      return std::make_pair(result, i + 1u);

    result *= ~schreier_structures[i]->transversal(beta);
  }

  return std::make_pair(result, base.size() + 1u);
}

bool BSGS::strips_completely(Perm const &perm) const
{
  auto strip_result(strip(perm));

  return strip_result.first.id() && strip_result.second == base.size() + 1u;
}

void BSGS::extend_base(unsigned bp)
{
  base.push_back(bp);
  schreier_structures.emplace_back(std::make_shared<SchreierTree>(degree()));
}

void BSGS::remove_generators()
{
  Dbg(Dbg::DBG) << "Removing redundant strong generators from BSGS:";
  Dbg(Dbg::DBG) << *this;

  Dbg(Dbg::TRACE) << "Stabilizers are:";
#ifndef NDEBUG
  for (auto i = 0u; i < base.size(); ++i)
    Dbg(Dbg::TRACE) << "S(" << i + 1u << ") = " << stabilizers(i);
#endif

  assert(base.size() > 0u);

  unsigned n = strong_generators[0].degree();
  unsigned k = base.size();

  std::unordered_set<Perm> strong_generator_set(
    strong_generators.begin(), strong_generators.end());


  auto difference = [&](std::unordered_set<Perm> const &lhs,
                        std::unordered_set<Perm> const &rhs) {

    int delta = static_cast<int>(lhs.size() - rhs.size());
    std::unordered_set<Perm> res;

    for (auto const &perm : lhs) {
      if (rhs.find(perm) == rhs.end() &&
          strong_generator_set.find(perm) != strong_generator_set.end()) {

        res.insert(perm);
        if (--delta == 0)
          break;
      }
    }

    return res;
  };

  auto produces_orbit = [&](unsigned root,
                            std::unordered_set<Perm> const &generators,
                            std::vector<unsigned> const &orbit) {

    std::vector<int> in_orbit_ref(n + 1u, 0), in_orbit(n + 1u, 0);

    for (unsigned x : orbit)
      in_orbit_ref[x] = 1;

    if (!in_orbit_ref[root])
      return false;

    in_orbit[root] = 1;

    std::vector<unsigned> queue {root};
    auto target = orbit.size() - 1u;

    while (!queue.empty()) {
      unsigned x = queue.back();
      queue.pop_back();

      for (auto const &gen : generators) {
        unsigned y = gen[x];
        if (!in_orbit_ref[y])
          return false;

        if (in_orbit[y] == 0) {
          in_orbit[y] = 1;
          queue.push_back(y);

          if (--target == 0u) {
            return true;
          }
        }
      }
    }

    return false;
  };

  std::unordered_set<Perm> stabilizer_set;
  std::unordered_set<Perm> stabilizer_intersection;

  for (int i = static_cast<int>(k - 1u); i >= 0; --i) {
    auto tmp(stabilizers(i));
    std::unordered_set<Perm> stabilizer_set_next(tmp.begin(), tmp.end());

    stabilizer_intersection = difference(stabilizer_set_next, stabilizer_set);

    stabilizer_set = stabilizer_set_next;

    Dbg(Dbg::TRACE) << "=== Considering S(" << i + 1u << ")/S(" << i + 2u << ")"
                    << " = " << stabilizer_intersection;

    if (stabilizer_intersection.size() < 2u)
      continue;

    auto it(stabilizer_intersection.begin());
    while (it != stabilizer_intersection.end()) {
      Dbg(Dbg::TRACE) << "Considering " << *it;

#ifndef NDEBUG
      std::unordered_set<Perm> reduced_stabilizers(stabilizer_set);

      auto it_(reduced_stabilizers.begin());
      while (it_ != reduced_stabilizers.end()) {
        if (*it_ == *it) {
          reduced_stabilizers.erase(it_++);
          break;
        } else {
          ++it_;
        }
      }
#endif
      auto orbit_gens(stabilizer_set);
      orbit_gens.erase(*it);

      bool remove_stab = false;
      if (produces_orbit(base[i], orbit_gens, orbit(i))) {
#ifndef NDEBUG
        Dbg(Dbg::TRACE) << base[i] << "^" << reduced_stabilizers
                        << " = " << orbit(i);
#endif

        Dbg(Dbg::TRACE) << "=> Removing strong generator " << *it;

        remove_stab = true;
      }
#ifndef NDEBUG
      else {
        Dbg(Dbg::TRACE) << base[i] << "^" << reduced_stabilizers
                        << " =/= " << orbit(i);
      }
#endif

      if (remove_stab) {
        strong_generator_set.erase(*it);
        stabilizer_set.erase(*it);
#ifndef NDEBUG
        reduced_stabilizers.erase(*it);
#endif
        stabilizer_intersection.erase(it++);
      } else {
        ++it;
      }
    }
  }

  strong_generators = std::vector<Perm>(
    strong_generator_set.begin(), strong_generator_set.end());

  Dbg(Dbg::TRACE) << "Remaining strong generators: " << strong_generators;

  std::vector<int> stabilizes(strong_generators.size(), 0);
  std::vector<Perm> stabilizers;

  for (int i = static_cast<int>(base.size()) - 1; i >= 0; --i) {
    for (auto j = 0u; j < strong_generators.size(); ++j) {
      if (stabilizes[j])
        continue;

      bool stab = true;
      for (int k = 0; k < i; ++k) {
        if (strong_generators[j][base[k]] != base[k]) {
          stab = false;
          break;
        }
      }

      if (stab) {
        stabilizes[j] = 1;
        stabilizers.push_back(strong_generators[j]);
      }
    }

    Dbg(Dbg::TRACE) << "=> S(" << i + 1u << ") = " << stabilizers;

    schreier_sims::orbit(base[i], stabilizers, schreier_structures[i]);
  }
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

    if (!bsgs.strips_completely(g)) {
      Dbg(Dbg::TRACE) << "=> Not in current BSGS";

      for (auto const &h : queue2) {
        Perm tmp(~g * ~h * g * h);
        if (!original_bsgs.strips_completely(tmp)) {
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

    while (!bsgs.strips_completely(gen)) {
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

std::ostream& operator<<(std::ostream& stream, BSGS const &bsgs)
{
  stream << "BASE: [";

  for (auto i = 0u; i < bsgs.base.size(); ++i) {
    stream << bsgs.base[i];
    if (i < bsgs.base.size() - 1u)
      stream << ", ";
  }

  stream << "]; SGS: [";

  for (auto i = 0u; i < bsgs.strong_generators.size(); ++i) {
    stream << bsgs.strong_generators[i];
    if (i < bsgs.strong_generators.size() - 1u)
      stream << ", ";
  }

  stream << ']';
  return stream;
}

} // namespace cgtl
