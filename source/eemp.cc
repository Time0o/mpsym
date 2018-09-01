#include <algorithm>
#include <climits>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>

#include "boost/container_hash/hash.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/strong_components.hpp"

#include "dbg.h"
#include "eemp.h"
#include "partial_perm.h"

namespace cgtl
{

std::vector<std::vector<unsigned>> EEMP::action_components(
  std::vector<unsigned> const &alpha, std::vector<PartialPerm> const &generators,
  SchreierTree &schreier_tree, OrbitGraph &orbit_graph)
{
#ifndef NDEBUG
  Dbg(Dbg::TRACE) << "Computing component of action of generators:";
  for (PartialPerm const &gen : generators)
    Dbg(Dbg::TRACE) << gen;

  Dbg(Dbg::TRACE) << "on:";
  Dbg(Dbg::TRACE) << alpha;
#endif

  // compute domain range over all generators
  unsigned dom_min = UINT_MAX;
  unsigned dom_max = 0u;
  for (PartialPerm const &gen : generators) {
    dom_min = std::min(dom_min, gen.dom_min());
    dom_max = std::max(dom_max, gen.dom_max());
  }
  unsigned dom_range_max = dom_max - dom_min + 1u;

  // data structures for efficient component membership testing
  std::vector<int> elem_size_present(dom_max + 1u, 0);
  elem_size_present[alpha.size()] = true;

  struct ComponentElemId {
    ComponentElemId(unsigned id, std::size_t hash) : id(id), hash(hash) {}

    unsigned id;
    std::size_t hash;
  };

  std::vector<std::vector<ComponentElemId>> elem_ids(dom_max + 1u);

  std::size_t alpha_hash = boost::hash_range(alpha.begin(), alpha.end());
  elem_ids[alpha.size()].push_back(ComponentElemId(0u, alpha_hash));

  // if beta is a component element, set 'id' to it's component index,
  // otherwise add add it to the component
  auto in_component = [&](std::vector<unsigned> const &beta, unsigned &id) {
    bool contained = true;

    auto const size_idx = beta.size();
    if (size_idx == dom_range_max && elem_size_present[size_idx]) {
      id = elem_ids[size_idx][0].id;
      return true;
    }

    std::size_t beta_hash = boost::hash_range(beta.begin(), beta.end());

    if (!elem_size_present[size_idx]) {
      contained = false;

    } else {
      auto it = std::find_if(
        elem_ids[size_idx].begin(), elem_ids[size_idx].end(),
        [&](ComponentElemId const &cei) { return cei.hash == beta_hash; });

      if (it == elem_ids[size_idx].end())
        contained = false;
      else
        id = it->id;
    }

    if (contained)
      return true;

    Dbg(Dbg::TRACE) << "Adjoining " << beta;
    elem_size_present[size_idx] = true;
    elem_ids[size_idx].push_back(ComponentElemId(id, beta_hash));

    return false;
  };


  // main loop
  std::vector<std::vector<unsigned>> component {alpha};
  decltype(schreier_tree.data) schreier_tree_data;
  decltype(orbit_graph.data) orbit_graph_data(generators.size());

  unsigned i = 0u, n = 0u;
  while (i < component.size()) {
    std::vector<unsigned> beta = component[i];
    Dbg(Dbg::TRACE) << "=== Looking at Component element " << beta
                    << " (i = " << i  << ')';

    for (auto j = 0u; j < generators.size(); ++j) {
      Dbg(Dbg::TRACE) << "== Considering generator " << generators[j];

      std::vector<unsigned> beta_prime = generators[j].image(beta);

      unsigned id = static_cast<unsigned>(component.size());

      if (!in_component(beta_prime, id)) {
        component.push_back(beta_prime);

        ++n;

        schreier_tree_data.push_back(std::make_pair(i, j));
        Dbg(Dbg::TRACE) << "v_" << n << " = " << i;
        Dbg(Dbg::TRACE) << "w_" << n << " = " << j;

        orbit_graph_data[j].push_back(n);
        Dbg(Dbg::TRACE) << "g_" << i << "," << j << " = " << n;

      } else {
        Dbg(Dbg::TRACE) << beta_prime << " already processed";

        orbit_graph_data[j].push_back(id);
        Dbg(Dbg::TRACE) << "g_" << i << "," << j << " = " << id;
      }
    }

    ++i;
  }

#ifndef NDEBUG
  std::stringstream ss;
  auto cellwidth = std::to_string(component.size()).length();

  // print resulting schreier tree
  ss << "i   |";
  ss << std::setw(cellwidth) << 2u;
  for (auto i = 1u; i < schreier_tree_data.size(); ++i)
    ss << ' ' << std::setw(cellwidth) << i + 2u;
  ss << '\n';

  ss << std::string(4u, '-');
  for (auto i = 0u; i < schreier_tree_data.size(); ++i)
    ss << std::string(cellwidth + 1u, '-');
  ss << '\n';

  ss << "v_i |";
  ss << std::setw(cellwidth) << std::get<0>(schreier_tree_data[0]) + 1u;
  for (auto i = 1u; i < schreier_tree_data.size(); ++i)
    ss << ' ' << std::setw(cellwidth) << std::get<0>(schreier_tree_data[i]) + 1u;
  ss << '\n';

  ss << "w_i |";
  ss << std::setw(cellwidth) << std::get<1>(schreier_tree_data[0]) + 1u;
  for (auto i = 1u; i < schreier_tree_data.size(); ++i)
    ss << ' ' << std::setw(cellwidth) << std::get<1>(schreier_tree_data[i]) + 1u;
  ss << '\n';

  Dbg(Dbg::TRACE) << "Resulting 'schreier tree':\n" << ss.str();

  // reset stringstream
  ss.str("");
  ss.clear();

  // print resulting orbit graph
  auto pad = std::to_string(orbit_graph_data.size()).length();

  ss << "i   " << std::string(pad, ' ') << " |";
  ss << std::setw(cellwidth) << 1u;
  for (auto i = 1u; i < orbit_graph_data[0].size(); ++i)
    ss << ' ' << std::setw(cellwidth) << i + 1u;
  ss << '\n';

  ss << std::string(5u + pad, '-');
  for (auto i = 0u; i < orbit_graph_data[0].size(); ++i)
    ss << std::string(cellwidth + 1, '-');
  ss << '\n';

  for (auto j = 0u; j < orbit_graph_data.size(); ++j) {
    ss << "g_i," << j << " |";

    ss << std::setw(cellwidth) << orbit_graph_data[j][0] + 1u;
    for (auto k = 1u; k < orbit_graph_data[j].size(); ++k)
      ss << ' ' << std::setw(cellwidth) << orbit_graph_data[j][k] + 1u;
    ss << '\n';
  }

  Dbg(Dbg::TRACE) << "Resulting orbit graph:\n" << ss.str();
#endif

  schreier_tree.dom_max = dom_max;
  schreier_tree.data = schreier_tree_data;

  orbit_graph.dom_max = dom_max;
  orbit_graph.data = orbit_graph_data;

  return component;
}

PartialPerm EEMP::schreier_trace(
  SchreierTree const &schreier_tree, unsigned i,
  std::vector<PartialPerm> const &generators)
{
  std::vector<unsigned> pperm(schreier_tree.dom_max);
  for (unsigned i = 0u; i < schreier_tree.dom_max; ++i)
    pperm[i] = i + 1u;

  PartialPerm res(pperm);

  while (i > 0u) {
    unsigned v = std::get<0>(schreier_tree.data[i - 1u]);
    unsigned w = std::get<1>(schreier_tree.data[i - 1u]);

    res = generators[w] * res;
    i = v;
  }

  return res;
}

std::vector<std::vector<unsigned>> EEMP::strongly_connected_components(
    OrbitGraph const &orbit_graph)
{
  boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> g;

  for (auto i = 0u; i < orbit_graph.data[0].size(); ++i)
    boost::add_vertex(g);

  for (auto i = 0u; i < orbit_graph.data[0].size(); ++i) {
    for (auto const &row : orbit_graph.data) {
      auto j = row[i];
      if (j != i)
        boost::add_edge(i, j, g);
    }
  }

  std::vector<int> component(boost::num_vertices(g));
  auto num = boost::strong_components(g, &component[0]);

  std::vector<std::vector<unsigned>> res(num);
  for (auto i = 0u; i < component.size(); ++i)
    res[component[i]].push_back(i);

  return res;
}

} // namespace cgtl
