#include <algorithm>
#include <functional>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/random_spanning_tree.hpp"

#include "dbg.h"
#include "eemp.h"
#include "partial_perm.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm.h"
#include "util.h"

namespace cgtl
{

PartialPermInverseSemigroup::PartialPermInverseSemigroup() : _empty(true) {}

PartialPermInverseSemigroup::PartialPermInverseSemigroup(
  std::vector<PartialPerm> const &generators) : _generators(generators)
{
  if (_generators.empty()) {
    _empty = true;
    return;
  } else {
    _empty = false;
  }

  unsigned dom_max = 0u;
  for (PartialPerm const &gen : generators)
    dom_max = std::max(dom_max, gen.dom_max());

  for (unsigned i = 0u; i < dom_max; ++i)
    _dom.push_back(i + 1u);

  std::vector<PartialPerm> inverse_generators(generators.size());
  for (auto i = 0u; i < generators.size(); ++i)
    inverse_generators[i] = ~generators[i];

  _ac_im = EEMP::action_component(_dom, generators, dom_max, _st_im, _og_im);

  EEMP::SchreierTree st_dummy;
  EEMP::OrbitGraph og_dummy;
  auto ac_dom(EEMP::action_component(
    _dom, inverse_generators, dom_max, st_dummy, og_dummy));

  for (auto i = 0u; i < _ac_im.size(); ++i)
    _ac_im_ht[_ac_im[i]] = i;

  _ac_dom_set = std::unordered_set<std::vector<unsigned>, VectorHash>(
    ac_dom.begin(), ac_dom.end());

  auto tmp(EEMP::strongly_connected_components(_og_im));
  unsigned num_scc = tmp.first;
  _scc = tmp.second;

  _scc_repr = std::vector<SccRepr>(num_scc);
  std::vector<int> found_repr(num_scc, 0);

  for (unsigned i = 0u; i < _scc.size(); ++i) {
    unsigned c = _scc[i];
    if (!found_repr[c]) {
      EEMP::SchreierTree st;
      EEMP::OrbitGraph og;
      auto ac(EEMP::action_component(
        _ac_im[i], _generators, _dom.back(), st, og));

      auto sg(EEMP::schreier_generators(
        _ac_im[i], _generators, _dom.back(), ac, st, og));

      _scc_repr[c] = SccRepr(i, _og_im, _scc, sg);
      found_repr[c] = 1;
    }
  }

  _r_class_repr = EEMP::r_class_representatives(_st_im, _generators);
}

bool PartialPermInverseSemigroup::is_element(PartialPerm const &pperm) const
{
  Dbg(Dbg::DBG) << "Testing membership of:";
  Dbg(Dbg::DBG) << pperm;
  Dbg(Dbg::DBG) << "in inverse semigroup with generators:";
  Dbg(Dbg::DBG) << _generators;

  if (_empty) {
    Dbg(Dbg::TRACE) << "Inverse semigroup is empty";
    Dbg(Dbg::DBG) << "==> Not an element";
    return false;
  }

  auto im(pperm.im());
  Dbg(Dbg::TRACE) << "Image is: " << im;

  auto ac_im_it(_ac_im_ht.find(im));
  if (ac_im_it == _ac_im_ht.end()) {
    Dbg(Dbg::TRACE) << "Image not compatible";
    Dbg(Dbg::DBG) << "==> Not an element";
    return false;
  }

  auto dom(pperm.dom());
  Dbg(Dbg::TRACE) << "Domain is: " << dom;

  auto it(std::find(_ac_dom_set.begin(), _ac_dom_set.end(), dom));
  if (it == _ac_dom_set.end()) {
    Dbg(Dbg::TRACE) << "Domain not compatible";
    Dbg(Dbg::DBG) << "==> Not an element";
    return false;
  }

  unsigned i = (*ac_im_it).second;
  SccRepr const &z_n = _scc_repr[_scc[i]];

  auto scc_repr(_ac_im[z_n.i]);
  Dbg(Dbg::TRACE) << "s.c.c representative is: " << scc_repr;

  auto vm_it(z_n.vertex_map.find(i));
  assert(vm_it != z_n.vertex_map.end());
  PartialPerm u(EEMP::schreier_trace((*vm_it).second, z_n.spanning_tree,
                                     _generators, _dom.back()));

  Dbg(Dbg::TRACE) << scc_repr << " * " << u << " = " << im;

  Dbg(Dbg::TRACE) << "=== Iterating over R class representatives:";
  Dbg(Dbg::TRACE) << "SGS of Sx is: " << z_n.schreier_generators.bsgs().sgs();
  for (PartialPerm const &x : _r_class_repr) {
    Dbg(Dbg::TRACE) << x;

    if (x.im() != scc_repr) {
      Dbg(Dbg::TRACE) << "=> Image not compatible";
      continue;
    }

    PartialPerm tmp_pperm(~x * pperm * ~u);
    if (tmp_pperm.id()) {
      Dbg(Dbg::TRACE) << "=> " << tmp_pperm << " is identity";
      Dbg(Dbg::DBG) << "==> element";
      return true;
    }

    Perm tmp_perm(tmp_pperm.to_perm(z_n.schreier_generators.degree()));

    if (z_n.schreier_generators.is_element(tmp_perm)) {
      Dbg(Dbg::TRACE) << "=> " << tmp_perm << " is contained in Sx";
      Dbg(Dbg::DBG) << "==> element";
      return true;
    }
#ifndef NDEBUG
    else
      Dbg(Dbg::TRACE) << "=> " << tmp_perm << " is not contained in Sx";
#endif
  }

  Dbg(Dbg::TRACE) << "Exhausted R class representatives";
  Dbg(Dbg::DBG) << "==> Not an element";
  return false;
}

PartialPermInverseSemigroup::SccRepr::SccRepr(unsigned i,
  EEMP::OrbitGraph const &orbit_graph, std::vector<unsigned> const &scc,
  PermGroup const &schreier_generators)
    : i(i), schreier_generators(schreier_generators)
{
  // construct s.c.c. subgraph
  boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> g;

  std::unordered_map<unsigned, unsigned> vertex_map_rev;

  for (unsigned j = 0u; j < scc.size(); ++j) {
    if (scc[j] == scc[i]) {
      unsigned v = static_cast<unsigned>(boost::add_vertex(g));
      vertex_map[j] = v;
      vertex_map_rev[v] = j;
    }
  }

  for (auto const &source : vertex_map) {
    for (auto const &row : orbit_graph.data) {
      auto dest(vertex_map.find(row[source.first]));
      if (dest != vertex_map.end() && (*dest).first != source.first)
        boost::add_edge((*dest).second, source.second, g);
    }
  }

  // find random spanning tree
  std::vector<unsigned> pred(boost::num_vertices(g));

  std::default_random_engine re;

  boost::random_spanning_tree(
    g, re, boost::root_vertex(vertex_map[i]).predecessor_map(&pred[0]));

  spanning_tree.data =
    std::vector<std::pair<unsigned, unsigned>>(pred.size() - 1u);

  for (auto j = 1u; j < pred.size(); ++j) {
    unsigned child = vertex_map_rev[j];
    unsigned parent = vertex_map_rev[pred[j]];

    // TODO: use edge properties instead
    for (unsigned gen = 0u; gen < orbit_graph.data.size(); ++gen) {
      if (orbit_graph.data[gen][parent] == child) {
        spanning_tree.data[j - 1u] = std::make_pair(pred[j], gen);
        break;
      }
    }
  }
}

} // namespace cgtl
