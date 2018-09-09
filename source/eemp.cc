#include <algorithm>
#include <climits>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <ostream>
#include <random>
#include <sstream>
#include <utility>
#include <vector>

#include "boost/container_hash/hash.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/strong_components.hpp"
#include "boost/graph/random_spanning_tree.hpp"

#include "dbg.h"
#include "eemp.h"
#include "partial_perm.h"
#include "perm.h"
#include "perm_group.h"
#include "util.h"

namespace cgtl
{

std::vector<std::vector<unsigned>> EEMP::action_component(
  std::vector<unsigned> const &alpha,
  std::vector<PartialPerm> const &generators, unsigned dom_max,
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
  for (PartialPerm const &gen : generators)
    dom_min = std::min(dom_min, gen.dom_min());

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

  std::size_t hash = boost::hash_range(alpha.begin(), alpha.end());
  elem_ids[alpha.size()].push_back(ComponentElemId(0u, hash));

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
        Dbg(Dbg::TRACE) << "=> v_" << n + 1u << " = " << i + 1u;
        Dbg(Dbg::TRACE) << "=> w_" << n  + 1u<< " = " << j + 1u;

        orbit_graph_data[j].push_back(n);
        Dbg(Dbg::TRACE) << "=> g_" << i + 1u << "," << j + 1u << " = " << n + 1u;

      } else {
        Dbg(Dbg::TRACE) << beta_prime << " already processed";

        orbit_graph_data[j].push_back(id);
        Dbg(Dbg::TRACE) << "=> g_" << i + 1u << "," << j + 1u << " = " << id + 1u;
      }
    }

    ++i;
  }

  schreier_tree.data = schreier_tree_data;

#ifndef NDEBUG
  Dbg(Dbg::TRACE) << "Resulting action component";
  for (auto i = 0u; i < component.size(); ++i)
    Dbg(Dbg::TRACE) << i + 1u << ": " << component[i];
#endif

  Dbg(Dbg::TRACE) << "Resulting schreier tree:\n" << schreier_tree;

  orbit_graph.data = orbit_graph_data;

  Dbg(Dbg::TRACE) << "Resulting orbit graph:\n" << orbit_graph;

  return component;
}

std::pair<unsigned, std::vector<unsigned>> EEMP::strongly_connected_components(
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

  std::vector<unsigned> component(boost::num_vertices(g));
  unsigned num = boost::strong_components(g, &component[0]);

  return std::make_pair(num, component);
}


EEMP::SchreierTree EEMP::scc_spanning_tree(
  unsigned i, OrbitGraph const &orbit_graph, std::vector<unsigned> const &scc)
{
  static std::default_random_engine re;

  Dbg(Dbg::TRACE) << "Finding spanning tree for s.c.c rooted at node " << i + 1u
                  << " in orbit graph:\n" << orbit_graph;

  // construct s.c.c. subgraph
  struct VertexProperty { unsigned i; };
  struct EdgeProperty { unsigned gen; };

  boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS, VertexProperty, EdgeProperty> g;

  std::vector<unsigned> scc_i;
  std::vector<unsigned> vertex_map(scc.size(), 0u);

  Dbg(Dbg::TRACE) << "== Constructing temporary subgraph";

  for (unsigned j = 0u; j < scc.size(); ++j) {
    if (scc[j] == scc[i]) {
      scc_i.push_back(j);

      unsigned v = boost::add_vertex({j}, g);
      vertex_map[j] = v + 1u;

      Dbg(Dbg::TRACE) << "Added vertex " << v << " (" << j + 1u << ')';
    }
  }

  assert(scc_i[0] == i && "provided representative is first node in s.c.c");

  for (unsigned source : scc_i) {
    for (unsigned row = 0; row < orbit_graph.data.size(); ++row) {
      unsigned dest = orbit_graph.data[row][source];

      unsigned source_vertex = vertex_map[source];
      unsigned dest_vertex = vertex_map[dest];

      if (dest_vertex > 0u && dest_vertex != source_vertex) {
        boost::add_edge(dest_vertex - 1u, source_vertex - 1u, {row}, g);

        Dbg(Dbg::TRACE) << "Added edge "
                        << source_vertex - 1u << " => " << dest_vertex - 1u
                        << " (" << source + 1u << " => " << dest + 1u
                        << " (" << row + 1u << "))";
      }
    }
  }

  // find random spanning tree
  SchreierTree spanning_tree;

  std::vector<unsigned> pred(boost::num_vertices(g));

  Dbg(Dbg::TRACE) << "Finding random spanning tree...";

  boost::random_spanning_tree(
    g, re, boost::root_vertex(vertex_map[i] - 1u).predecessor_map(&pred[0]));

  Dbg(Dbg::TRACE) << "Found spanning tree with predecessor relationships: "
                  << std::vector<unsigned>(pred.begin() + 1, pred.end());

  // convert into Schreier tree format
  spanning_tree.data =
    std::vector<std::pair<unsigned, unsigned>>(scc.size() - 1u);

  for (auto j = 1u; j < pred.size(); ++j) {
    auto child(boost::vertex(j, g));
    auto parent(boost::vertex(pred[j], g));
    auto edge(boost::edge(child, parent, g));

    assert(edge.second && "edge exists");

    spanning_tree.data[g[child].i - 1u] =
      std::make_pair(g[parent].i, g[edge.first].gen);
  }

  Dbg(Dbg::TRACE) << "Resulting spanning Schreier tree is:\n" << spanning_tree;

  return spanning_tree;
}

PartialPerm EEMP::schreier_trace(
  unsigned x, SchreierTree const &schreier_tree,
  std::vector<PartialPerm> const &generators, unsigned dom_max, unsigned target)
{
  PartialPerm res(dom_max);

  while (x != target) {
    unsigned v = std::get<0>(schreier_tree.data[x - 1u]);
    unsigned w = std::get<1>(schreier_tree.data[x - 1u]);

    res = generators[w] * res;
    x = v;
  }

  return res;
}

PermGroup EEMP::schreier_generators(unsigned x,
  std::vector<PartialPerm> const &generators, unsigned dom_max,
  std::vector<std::vector<unsigned>> const &action_component,
  SchreierTree const &schreier_tree, OrbitGraph const &orbit_graph,
  std::vector<unsigned> const &sccs)
{
  auto im(action_component[x]);

  Dbg(Dbg::TRACE) << "Finding schreier generators for Sx for: " << im;

  if (im.empty()) {
    Dbg(Dbg::TRACE) << "==> Returning empty permutation group";
    return PermGroup();
  }

  unsigned im_max = im.back();

  Dbg(Dbg::TRACE) << "Finding schreier generators for Sx for: " << im;

  std::vector<unsigned> scc;
  for (unsigned i = 0u; i < sccs.size(); ++i) {
    if (sccs[i] == sccs[x])
      scc.push_back(i);
  }

#ifndef NDEBUG
  std::vector<std::vector<unsigned>> _scc(scc.size());
  for (auto i = 0u; i < scc.size(); ++i)
    _scc[i] = action_component[scc[i]];

  Dbg(Dbg::TRACE) << "Strongly connected component of x is: " << _scc;
#endif

  std::vector<Perm> sg_gens;

  for (auto i = 0u; i < scc.size(); ++i) {
    for (auto j = 0u; j < generators.size(); ++j) {

      unsigned k = orbit_graph.data[j][scc[i]];
      if (sccs[k] != sccs[x])
        continue;

      Dbg(Dbg::TRACE) << "Found generator for j/i(j)/k/l = "
                      << i + 1u << '/'
                      << scc[i] + 1u << '/'
                      << j + 1u << '/'
                      << k + 1u;

      PartialPerm ui(schreier_trace(
        scc[i], schreier_tree, generators, dom_max, x));
      PartialPerm uk(schreier_trace(
        k, schreier_tree, generators, dom_max, x));

      Dbg(Dbg::TRACE) << "where:";
      Dbg(Dbg::TRACE) << "u_i(j) = " << "u_" << scc[i] + 1u << " = " << ui;
      Dbg(Dbg::TRACE) << "x_k = " << "x_" << j + 1u << " = " << generators[j];
      Dbg(Dbg::TRACE) << "~u_l = " << "u_" << k + 1u << " = " << ~uk;

      PartialPerm sg(ui * generators[j] * ~uk);
      sg = sg.restricted(im);

      Dbg(Dbg::TRACE) << "=> Schreier generator is: " << sg;

      if (!sg.id())
        sg_gens.push_back(sg.to_perm(im_max));
    }
  }

  PermGroup res(im_max, sg_gens);

  Dbg(Dbg::TRACE) << "==> Returning:";
  Dbg(Dbg::TRACE) << res;

  return res;
}

std::vector<PartialPerm> EEMP::r_class_representatives(
  SchreierTree const &schreier_tree, std::vector<PartialPerm> const &generators)
{
  std::vector<std::vector<unsigned>> st_adj(schreier_tree.data.size() + 1u);
  for (auto i = 0u; i < schreier_tree.data.size(); ++i)
    st_adj[std::get<0>(schreier_tree.data[i])].push_back(i + 1u);

  std::vector<PartialPerm> res;

  std::function<void(unsigned, PartialPerm const &)>
  backtrack = [&](unsigned node, PartialPerm const &pperm) {
    res.push_back(pperm);

    for (unsigned child : st_adj[node]) {
      unsigned gen_idx = schreier_tree.data[child - 1u].second;
      PartialPerm next_pperm(pperm * generators[gen_idx]);

      backtrack(child, next_pperm);
    }
  };

  backtrack(0, generators[0]);

  return res;
}

std::ostream& operator<<(
  std::ostream& stream, EEMP::SchreierTree const &schreier_tree)
{
  if (schreier_tree.data.empty()) {
    stream << "empty schreier tree";
    return stream;
  }

  auto size = schreier_tree.data.size();
  auto cellwidth = std::to_string(size).length();

  stream << "i   |";
  stream << std::setw(cellwidth) << 2u;
  for (auto i = 1u; i < size; ++i)
    stream << ' ' << std::setw(cellwidth) << i + 2u;
  stream << '\n';

  stream << std::string(4u, '-');
  for (auto i = 0u; i < size; ++i)
    stream << std::string(cellwidth + 1u, '-');
  stream << '\n';

  stream << "v_i |";
  stream << std::setw(cellwidth) << std::get<0>(schreier_tree.data[0]) + 1u;
  for (auto i = 1u; i < size; ++i) {
    stream << ' ' << std::setw(cellwidth)
           << std::get<0>(schreier_tree.data[i]) + 1u;
  }
  stream << '\n';

  stream << "w_i |";
  stream << std::setw(cellwidth) << std::get<1>(schreier_tree.data[0]) + 1u;
  for (auto i = 1u; i < size; ++i) {
    stream << ' ' << std::setw(cellwidth)
           << std::get<1>(schreier_tree.data[i]) + 1u;
  }

  return stream;
}

std::ostream& operator<<(
  std::ostream& stream, EEMP::OrbitGraph const &orbit_graph)
{
  if (orbit_graph.data.empty()) {
    stream << "empty orbit graph";
    return stream;
  }

  auto size = orbit_graph.data[0].size();

  auto cellwidth = std::to_string(size).length();
  auto pad = std::to_string(orbit_graph.data.size()).length();

  stream << "i   " << std::string(pad, ' ') << "  |";
  stream << std::setw(cellwidth) << 1u;
  for (auto i = 1u; i < size; ++i)
    stream << ' ' << std::setw(cellwidth) << i + 1u;
  stream << '\n';

  stream << std::string(6u + pad, '-');
  for (auto i = 0u; i < size; ++i)
    stream << std::string(cellwidth + 1, '-');
  stream << '\n';

  for (auto j = 0u; j < orbit_graph.data.size(); ++j) {
    stream << "g_i," << std::setw(pad) << j + 1u << "  |";

    stream << std::setw(cellwidth) << orbit_graph.data[j][0] + 1u;
    for (auto k = 1u; k < size; ++k)
      stream << ' ' << std::setw(cellwidth) << orbit_graph.data[j][k] + 1u;
    stream << '\n';
  }

  auto tmp(EEMP::strongly_connected_components(orbit_graph));
  auto scc(expand_partition(tmp.second));

  stream << "s.c.c." << std::string(pad - 1u, ' ') << " | {";
  for (auto i = 0u; i < scc.size(); ++i) {
    stream << '{' << scc[i][0] + 1u;
    for (auto j = 1u; j < scc[i].size(); ++j)
      stream << ", " << scc[i][j] + 1u;
    stream << '}';
    if (i != scc.size() - 1u)
      stream << ", ";
  }
  stream << '}';

  return stream;
}

} // namespace cgtl
