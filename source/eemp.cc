#include <algorithm>
#include <climits>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <utility>
#include <utility>
#include <vector>

#include "boost/container_hash/hash.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/strong_components.hpp"

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

std::vector<unsigned> EEMP::strongly_connected_components(
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
  boost::strong_components(g, &component[0]);
  return component;
}

PartialPerm EEMP::schreier_trace(
  unsigned x, SchreierTree const &schreier_tree,
  std::vector<PartialPerm> const &generators, unsigned dom_max)
{
  std::vector<unsigned> pperm(dom_max);
  for (unsigned j = 0u; j < dom_max; ++j)
    pperm[j] = j + 1u;

  PartialPerm res(pperm);

  while (x > 0u) {
    unsigned v = std::get<0>(schreier_tree.data[x - 1u]);
    unsigned w = std::get<1>(schreier_tree.data[x - 1u]);

    res = generators[w] * res;
    x = v;
  }

  return res;
}

PermGroup EEMP::schreier_generators(PartialPerm const &x,
  std::vector<PartialPerm> const &generators, unsigned dom_max,
  std::vector<std::vector<unsigned>> const &action_component,
  SchreierTree const &schreier_tree, OrbitGraph const &orbit_graph)
{
  Dbg(Dbg::TRACE) << "Finding schreier generators for Sx where x is: " << x;

  auto im(x.im());
  unsigned im_max = *std::max_element(im.begin(), im.end());

  if (im.empty())
    return PermGroup();

  auto sccs(strongly_connected_components(orbit_graph));

  std::vector<unsigned> scc;
  for (auto i = 0u; i < sccs.size(); ++i) {
    if (sccs[i] == sccs[0])
      scc.push_back(static_cast<unsigned>(i));
  }

#ifndef NDEBUG
  std::vector<std::vector<unsigned>> _scc(scc.size());
  for (auto i = 0u; i < scc.size(); ++i)
    _scc[i] = action_component[scc[i]];

  Dbg(Dbg::TRACE) << "Strongly connected component of x is: " << _scc;
#endif

  std::vector<Perm> res;

  for (auto i = 0u; i < scc.size(); ++i) {
    for (auto j = 0u; j < generators.size(); ++j) {

      unsigned k = orbit_graph.data[j][scc[i]];
      if (sccs[k] != sccs[0])
        continue;

      Dbg(Dbg::TRACE) << "Found generator for j/i(j)/k/l = "
                      << i + 1u << '/'
                      << scc[i] + 1u << '/'
                      << j + 1u << '/'
                      << k + 1u;

      PartialPerm ui(schreier_trace(scc[i], schreier_tree, generators, dom_max));
      PartialPerm uk(schreier_trace(k, schreier_tree, generators, dom_max));

      Dbg(Dbg::TRACE) << "where:";
      Dbg(Dbg::TRACE) << "u_i(j) = " << "u_" << scc[i] + 1u << " = " << ui;
      Dbg(Dbg::TRACE) << "x_k = " << "x_" << j + 1u << " = " << generators[j];
      Dbg(Dbg::TRACE) << "~u_l = " << "u_" << k + 1u << " = " << ~uk;

      PartialPerm sg(ui * generators[j] * ~uk);
      sg = sg.restricted(im);

      Dbg(Dbg::TRACE) << "=> Schreier generator is: " << sg;

      if (!sg.id())
        res.push_back(sg.to_perm(im_max));
    }
  }

  return PermGroup(im_max, res);
}

std::vector<PartialPerm> EEMP::r_class_elements(PartialPerm const &x,
  std::vector<PartialPerm> const &generators, unsigned dom_max)
{
  SchreierTree st;
  OrbitGraph og;
  auto ac(action_component(x.im(), generators, dom_max, st, og));
  auto sccs(strongly_connected_components(og));
  auto sx(schreier_generators(x, generators, dom_max, ac, st, og));

  std::vector<unsigned> scc;
  for (unsigned i = 0u; i < sccs.size(); ++i) {
    if (sccs[i] == sccs[0])
      scc.push_back(i);
  }

  std::vector<PartialPerm> res;

  for (unsigned i : scc) {
    for (Perm const &s : sx) {
      PartialPerm ui(schreier_trace(i, st, generators, dom_max));
      res.push_back(x * PartialPerm::from_perm(s) * ui);
    }
  }

  return res;
}

std::vector<PartialPerm> EEMP::r_classes_in_d_class(PartialPerm const &x,
  std::vector<PartialPerm> const &generators, unsigned dom_max)
{
  std::vector<PartialPerm> inverted_generators(generators.size());
  for (auto i = 0u; i < generators.size(); ++i)
    inverted_generators[i] =  ~generators[i];

  SchreierTree st_dom;
  OrbitGraph og_dom;
  auto ac_dom(action_component(
    x.dom(), inverted_generators, dom_max, st_dom, og_dom));

  auto xs(schreier_generators(
    ~x, inverted_generators, dom_max, ac_dom, st_dom, og_dom));

  auto scc(strongly_connected_components(og_dom));

  std::vector<std::vector<unsigned>> st_dom_adj(st_dom.data.size() + 1u);

  for (auto i = 0u; i < st_dom.data.size(); ++i)
    st_dom_adj[std::get<0>(st_dom.data[i])].push_back(i + 1u);

  std::vector<PartialPerm> res;

  std::function<void(unsigned, PartialPerm const &)>
  backtrack = [&](unsigned node, PartialPerm const &pperm) {
    res.push_back(pperm);

    for (unsigned child : st_dom_adj[node]) {
      unsigned gen_idx = std::get<1>(st_dom.data[child - 1u]);
      PartialPerm next_pperm(generators[gen_idx] * pperm);

      if (scc[child] == scc[0])
        backtrack(child, next_pperm);
    }
  };

  backtrack(0, x);

  return res;
}

std::vector<PartialPerm> EEMP::inverse_semigroup_r_class_representatives(
  EEMP::S const &s)
{
  std::vector<std::vector<unsigned>> st_adj(s.schreier_tree.data.size() + 1u);
  for (auto i = 0u; i < s.schreier_tree.data.size(); ++i)
    st_adj[std::get<0>(s.schreier_tree.data[i])].push_back(i + 1u);

  std::vector<PartialPerm> res;

  std::function<void(unsigned, PartialPerm const &)>
  backtrack = [&](unsigned node, PartialPerm const &pperm) {
    res.push_back(pperm);

    for (unsigned child : st_adj[node]) {
      unsigned gen_idx = std::get<1>(s.schreier_tree.data[child - 1u]);
      PartialPerm next_pperm(s.generators[gen_idx] * pperm);

      backtrack(child, next_pperm);
    }
  };

  backtrack(0, s.alpha);

  return res;
}

bool EEMP::is_inverse_semigroup_member(PartialPerm const &y,
  EEMP::S const &s, PermGroup const &sx)
{
  std::vector<PartialPerm> inverted_generators(s.generators.size());
  for (auto i = 0u; i < s.generators.size(); ++i)
    inverted_generators[i] = ~s.generators[i];

  SchreierTree left_schreier_tree;
  OrbitGraph left_orbit_graph;
  auto left_action_component(
    action_component(s.alpha, inverted_generators, s.dom,
                     left_schreier_tree, left_orbit_graph));

  auto y_dom(y.dom());
  auto y_im(y.im());

  auto right_it = std::find(s.action_component.begin(),
                            s.action_component.end(), y_im);

  auto left_it = std::find(left_action_component.begin(),
                           left_action_component.end(), y_dom);

  if (right_it == s.action_component.end() ||
      left_it == left_action_component.end()) {
    return false;
  }

  // find n / alpha_n
  std::vector<unsigned> alpha_n;

  int _n = -1;
  for (int i = 0u; i < static_cast<int>(s.scc.size()); ++i) {
    if (s.action_component[i] == y_im) {
      alpha_n = s.action_component[i];
      _n = s.scc[i];
      break;
    }
  }

  assert(_n >= 0);

  unsigned n = static_cast<unsigned>(_n);

  // find u
  PartialPerm u;

  unsigned y_idx = static_cast<unsigned>(right_it - s.action_component.begin());
  for (auto i = 0u; i < s.orbit_graph.data.size(); ++i) {
    if (s.orbit_graph.data[i][n] == y_idx) {
      u = s.generators[i];
      break;
    }
  }

  assert(!u.empty());

  for (PartialPerm const &x : inverse_semigroup_r_class_representatives(s)) {
     if (x.im() != alpha_n)
       continue;

     if (x.dom() != y.dom())
       continue;

     if (sx.is_element(PartialPerm(~x * y * ~u).to_perm(sx.degree())))
       return true;
  }

  return false;
}

EEMP::S::S(std::vector<PartialPerm> const &generators) : generators(generators)
{
  dom = 0u;
  for (auto const &gen : generators)
    dom = std::max(dom, gen.dom_max());

  std::vector<unsigned> alpha(dom);
  for (unsigned i = 0u; i < dom; ++i)
    alpha[i] = i + 1u;

  action_component = EEMP::action_component(
    alpha, generators, dom, schreier_tree, orbit_graph);

  scc = strongly_connected_components(orbit_graph);
}

EEMP::RClass::RClass(PartialPerm const &x, EEMP::S const &s)
  : _x(x), _s(s)
{
  _action_component = EEMP::action_component(_x.im(), _s.generators, _s.dom,
    _schreier_tree, _orbit_graph);

  _schreier_generators = EEMP::schreier_generators(_x, _s.generators, _s.dom,
    _action_component, _schreier_tree, _orbit_graph);
}

bool EEMP::RClass::is_member(PartialPerm const &y)
{
  Dbg(Dbg::TRACE) << "Testing R class membership:";
  Dbg(Dbg::TRACE) << "y = " << y;
  Dbg(Dbg::TRACE) << "x = " << _x;
  Dbg(Dbg::TRACE) << "S ~ " << _s.generators;

  if (y.dom() != _x.dom()) {
    Dbg(Dbg::TRACE) << "x and y have different domains => no member";
    return false;
  }

  auto x_im(_x.im());
  auto y_im(y.im());

  int ac_idx_x = -1;
  int ac_idx_y = -1;
  for (auto i = 0u; i < _s.action_component.size(); ++i) {
    if (ac_idx_x < 0 && _s.action_component[i] == x_im)
      ac_idx_x = i;

    if (ac_idx_y < 0 && _s.action_component[i] == y_im)
      ac_idx_y = i;

    if (ac_idx_x >= 0 && ac_idx_y >= 0)
      break;
  }

  assert(ac_idx_x >= 0 && ac_idx_y >= 0);

  if (_s.scc[ac_idx_x] != _s.scc[ac_idx_y]) {
    Dbg(Dbg::TRACE) << "x and y do not lie in same s.c.c => no member";
    return false;
  }

  PartialPerm u;
  for (auto i = 0u; i < _s.orbit_graph.data.size(); ++i) {
    if (_s.orbit_graph.data[i][ac_idx_x] == static_cast<unsigned>(ac_idx_y))
       u = _s.generators[i];
  }

  assert(!u.empty());

  Dbg(Dbg::TRACE) << "u = " << u;

  unsigned x_im_max = *std::max_element(x_im.begin(), x_im.end());
  Perm perm(PartialPerm(~_x * y * ~u).to_perm(x_im_max));

  bool res = _schreier_generators.is_element(perm);

  if (res)
    Dbg(Dbg::TRACE) << perm << "is a member of Sx => member";
  else
    Dbg(Dbg::TRACE) << perm << "is not a member of Sx => no member";

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

  auto scc(EEMP::strongly_connected_components(orbit_graph));
  auto scc_expanded(expand_partition(scc));

  stream << "s.c.c." << std::string(pad - 1u, ' ') << " | {";
  for (auto i = 0u; i < scc_expanded.size(); ++i) {
    stream << '{' << scc_expanded[i][0] + 1u;
    for (auto j = 1u; j < scc_expanded[i].size(); ++j)
      stream << ", " << scc_expanded[i][j] + 1u;
    stream << '}';
    if (i != scc_expanded.size() - 1u)
      stream << ", ";
  }
  stream << '}';

  return stream;
}

} // namespace cgtl
