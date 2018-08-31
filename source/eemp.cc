#include "dbg.h"
#include "eemp.h"
#include "partial_perm.h"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <vector>

#include "boost/container_hash/hash.hpp"

namespace cgtl
{

void EEMP::action_components(
  std::vector<unsigned> const &alpha, std::vector<PartialPerm> const &generators)
{
  Dbg(Dbg::TRACE) << "Computing component of action of " << generators
                  << " on "<< alpha;

  unsigned dom_max = 0u;
  for (PartialPerm const &gen : generators)
    dom_max = std::max(dom_max, gen.dom_max());

  std::vector<std::vector<unsigned>> component;

  std::vector<int> elem_size_present(dom_max, 0);
  std::vector<std::vector<std::size_t>> elem_hashes(dom_max);

  auto in_component = [&](std::vector<unsigned> const &beta) {
    bool contained = true;

    auto const size_idx = beta.size() - 1u;
    std::size_t beta_hash = boost::hash_range(beta.begin(), beta.end());

    if (!elem_size_present[size_idx]) {
      contained = false;

    } else {
      auto it = std::find(elem_hashes[size_idx].begin(),
                          elem_hashes[size_idx].end(), beta_hash);

      if (it == elem_hashes[size_idx].end())
        contained = false;
    }

    if (contained)
      return true;

    Dbg(Dbg::TRACE) << "Adjoining " << beta;
    elem_hashes[size_idx].push_back(beta_hash);

    return false;
  };

  component.push_back(alpha);

  struct Schreier { unsigned v, w; };
  std::vector<Schreier> schreier_tree(dom_max);

  std::vector<unsigned> orbit_graph(generators.size() * dom_max);

  unsigned i = 0u, n = 0u;
  while (i < component.size()) {
    std::vector<unsigned> beta = component[i];
    Dbg(Dbg::TRACE) << "=== Looking at Component element " << beta
                    << " (i = " << i  << ')';

    for (auto j = 0u; j < generators.size(); ++j) {
      PartialPerm gen = generators[i];
      Dbg(Dbg::TRACE) << "== Considering generator " << gen;

      std::vector<unsigned> beta_prime = gen.image(beta);

      if (!in_component(beta_prime)) {
        ++n;
        schreier_tree[n] = {i, j};
        orbit_graph[j * dom_max + i] = n;
      } else {
        orbit_graph[j * dom_max + i] = 0; // r???
      }
    }

    ++i;
  }

  Dbg(Dbg::TRACE) << "Resulting 'schreier tree':";

#ifndef NDEBUG
  std::stringstream ss;

  auto cellwidth = std::to_string(dom_max).length();
  ss << std::setw(cellwidth);

  ss << "i   |";
  for (unsigned i = 0u; i < dom_max; ++i)
    ss << i;
  ss << '\n';

  for (unsigned i = 0u; i < dom_max; ++i)
    ss << std::string(cellwidth, '-');
  ss << '\n';

  ss << "v_i |";
  for (unsigned i = 0u; i < dom_max; ++i)
    ss << schreier_tree[i].v;
  ss << '\n';

  ss << "w_i |";
  for (unsigned i = 0u; i < dom_max; ++i)
    ss << schreier_tree[i].w;
  ss << '\n';
#endif
}

} // namespace cgtl
