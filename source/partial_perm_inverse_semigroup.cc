#include <algorithm>
#include <cassert>
#include <functional>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random_spanning_tree.hpp>

#include "dbg.h"
#include "dump.h"
#include "eemp.h"
#include "partial_perm.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm.h"
#include "util.h"

/**
 * @file partial_perm_inverse_semigroup.cc
 * @author Timo Nicolai
 *
 * @brief Implements `PartialPermInverseSemigroup`.
 */

namespace mpsym
{

PartialPermInverseSemigroup::PartialPermInverseSemigroup() : _trivial(true) {}

PartialPermInverseSemigroup::PartialPermInverseSemigroup(
  std::vector<PartialPerm> const &generators)
    : _trivial(generators.empty()), _generators(generators)
{
  unsigned dom_max = 0u;
  for (PartialPerm const &gen : generators)
    dom_max = std::max(dom_max, gen.dom_max());

  for (unsigned i = 0u; i < dom_max; ++i)
    _dom.push_back(i + 1u);

  std::vector<PartialPerm> inverse_generators(generators.size());
  for (auto i = 0u; i < generators.size(); ++i)
    inverse_generators[i] = ~generators[i];

  _ac_im = eemp::action_component(_dom, generators, dom_max, _st_im, _og_im);

  eemp::SchreierTree st_dummy;
  eemp::OrbitGraph og_dummy;
  auto ac_dom(eemp::action_component(
    _dom, inverse_generators, dom_max, st_dummy, og_dummy));

  for (auto i = 0u; i < _ac_im.size(); ++i)
    _ac_im_ht[_ac_im[i]] = i;

  update_scc_representatives();

  _r_class_repr = eemp::r_class_representatives(_st_im, _generators);
}

bool PartialPermInverseSemigroup::contains_element(
  PartialPerm const &pperm) const
{
  DBG(DEBUG) << "Testing membership of:";
  DBG(DEBUG) << pperm;
  DBG(DEBUG) << "in inverse semigroup with generators:";
  DBG(DEBUG) << _generators;

  if (_trivial) {
    DBG(TRACE) << "Inverse semigroup is empty";
    if (pperm.empty()) {
      DBG(DEBUG) << "==> Not an element";
      return true;
    } else {
      DBG(DEBUG) << "==> Not an element";
      return false;
    }
  }

  auto im(pperm.im());
  DBG(TRACE) << "Image is: " << im;

  auto ac_im_it(_ac_im_ht.find(im));
  if (ac_im_it == _ac_im_ht.end()) {
    DBG(TRACE) << "Image not compatible";
    DBG(DEBUG) << "==> Not an element";
    return false;
  }

  auto dom(pperm.dom());
  DBG(TRACE) << "Domain is: " << dom;

  auto ac_dom_it(_ac_im_ht.find(dom));
  if (ac_dom_it == _ac_im_ht.end()) {
    DBG(TRACE) << "Domain not compatible";
    DBG(DEBUG) << "==> Not an element";
    return false;
  }

  unsigned i = (*ac_im_it).second;
  SccRepr const &z_n = _scc_repr[_scc[i]];

  auto scc_repr(_ac_im[z_n.i]);
  DBG(TRACE) << "s.c.c representative is: " << scc_repr;

  PartialPerm u(eemp::schreier_trace(i, z_n.spanning_tree,
                                     _generators, _dom.back(), z_n.i));

  DBG(TRACE) << scc_repr << " * " << u << " = " << im;

  DBG(TRACE) << "=== Iterating over R class representatives:";
  DBG(TRACE) << "SGS of Sx is: " << z_n.schreier_generators.generators();

  for (PartialPerm const &x : _r_class_repr) {
    DBG(TRACE) << x;

    if (x.im() != scc_repr) {
      DBG(TRACE) << "=> Image not compatible";
      continue;
    }

    PartialPerm tmp_pperm(~x * pperm * ~u);
    if (tmp_pperm.id()) {
      DBG(TRACE) << "=> " << tmp_pperm << " is identity";
      DBG(DEBUG) << "==> Element";
      return true;
    }

    Perm tmp_perm(tmp_pperm.to_perm(z_n.schreier_generators.degree()));

    if (z_n.schreier_generators.contains_element(tmp_perm)) {
      DBG(TRACE) << "=> " << tmp_perm << " is contained in Sx";
      DBG(DEBUG) << "==> Element";
      return true;
    }
#ifndef NDEBUG
    else
      DBG(TRACE) << "=> " << tmp_perm << " is not contained in Sx";
#endif
  }

  DBG(TRACE) << "Exhausted R class representatives";
  DBG(DEBUG) << "==> Not an element";
  return false;
}

void PartialPermInverseSemigroup::adjoin_generators(
  std::vector<PartialPerm> const &generators,
  bool minimize)
{
  DBG(DEBUG) << "Adjoining generators: " << generators;
  DBG(DEBUG) << "to inverse semigroups with generators: " << _generators;

  if (generators.empty())
    return;

  if (_trivial) {
    DBG(TRACE) << "==> Creating new inverse semigroup";
    *this = PartialPermInverseSemigroup(generators);
    return;
  }

  DBG(TRACE) << "=== Adjoining new generators";

  if (!minimize) {
    update_action_component(generators);

    _generators.insert(_generators.end(), generators.begin(), generators.end());
    DBG(TRACE) << "New generator set is: " << _generators;

  } else {
    DBG(TRACE) << "Trying to minimize number of generators adjoined";

    for (auto j = 0u; j < generators.size(); ++j) {
      auto const &gen(generators[j]);

      if (contains_element(gen)) {
        DBG(TRACE) << "Skipping " << gen << " (is already an element)";
        continue;
      }
      DBG(TRACE) << "Adjoining " << gen;

      if (gen.dom_max() > _dom.back()) {
        for (auto i = _dom.size() + 1u; i <= gen.dom_max(); ++i)
          _dom.push_back(i);
        DBG(TRACE) << "Extended domain to: " << _dom;
      }

      update_action_component({gen});

      _generators.push_back(gen);
      DBG(TRACE) << "New generator set is: " << _generators;
    }
  }

  DBG(TRACE) << "Updating s.c.c representatives";
  update_scc_representatives();

  DBG(TRACE) << "Updating R class representatives";
  _r_class_repr = eemp::r_class_representatives(_st_im, _generators);
}

void PartialPermInverseSemigroup::update_action_component(
  std::vector<PartialPerm> const &generators)
{
  DBG(TRACE) << "Initial action component:\n" << _ac_im;
  DBG(TRACE) << "Initial orbit graph:\n" << _og_im;
  DBG(TRACE) << "Initial Schreier tree:\n" << _st_im;

  unsigned first_new_row = _og_im.data.size();
  unsigned first_new_node = _og_im.data[0].size();
  unsigned next_new_node = first_new_node;

  for (auto i = 0u; i < generators.size(); ++i)
    _og_im.data.push_back(std::vector<unsigned>(_og_im.data[0].size()));

  std::queue<unsigned> node_queue;
  node_queue.push(0u);

  std::vector<int> visited(_og_im.data[0].size(), 0);

  DBG(TRACE) << "=== Extending orbit graph";

  while (!node_queue.empty()) {
    unsigned node = node_queue.front();
    node_queue.pop();

    DBG(TRACE) << "== Considering node: " << _ac_im[node];

    // add new column to orbit graph
    if (node >= _og_im.data[0].size()) {
      DBG(TRACE) << "Adding new column to orbit graph";
      for (auto row = 0u; row < _og_im.data.size(); ++row)
        _og_im.data[row].push_back(0u);
    }

    auto n_gens = generators.size();
    if (node >= first_new_node) {
      // also apply old generators to newly created nodes
      n_gens += _generators.size();
    } else {
      // otherwise make sure to add ALL children to the queue
      for (auto i = 0u; i < _generators.size(); ++i) {
        unsigned child = _og_im.data[i][node];
        if (!visited[child])
          node_queue.push(child);
      }
    }

    for (auto i = 0u; i < n_gens; ++i) {
      PartialPerm gen;
      unsigned row_idx;

      if (i >= generators.size()) {
        row_idx = i - generators.size();
        gen = _generators[row_idx];
      } else {
        row_idx = first_new_row + i;
        gen = generators[i];
      }

      DBG(TRACE) << "Applying generator: " << gen;

      auto next_node(gen.image<std::vector>(_ac_im[node].begin(),
                                            _ac_im[node].end()));

      DBG(TRACE) << "=> Next node is: " << next_node;

      auto it(_ac_im_ht.find(next_node));
      if (it == _ac_im_ht.end()) {
        DBG(TRACE) << "==> Adding new node/edge to orbit graph";
        _og_im.data[row_idx][node] = next_new_node;
        _ac_im.push_back(next_node);
        _ac_im_ht[next_node] = next_new_node;
        node_queue.push(next_new_node++);
        visited.push_back(0);

        DBG(TRACE) << "==> Updating Schreier tree";
        _st_im.data.push_back(std::make_pair(node, row_idx));
      } else {
        DBG(TRACE) << "==> Adding new edge to orbit graph";
        _og_im.data[row_idx][node] = (*it).second;
      }
    }

    visited[node] = 1;
  }

  DBG(TRACE) << "Updating s.c.c representatives";

  DBG(TRACE) << "Resulting action component:";
#ifndef NDEBUG
  for (auto const &c : _ac_im)
    DBG(TRACE) << c;
#endif
  DBG(TRACE) << "Resulting orbit graph:\n" << _og_im;
  DBG(TRACE) << "Resulting Schreier tree:\n" << _st_im;
}

void PartialPermInverseSemigroup::update_scc_representatives()
{
  auto tmp(eemp::strongly_connected_components(_og_im));
  unsigned num_scc = tmp.first;
  _scc = tmp.second;

  _scc_repr = std::vector<SccRepr>(num_scc);
  std::vector<int> found_repr(num_scc, 0);

  for (unsigned i = 0u; i < _scc.size(); ++i) {
    unsigned c = _scc[i];
    if (!found_repr[c]) {
      auto st(eemp::scc_spanning_tree(i, _og_im, _scc));

      auto sg(eemp::schreier_generators(
        i, _generators, _dom.back(), _ac_im, st, _og_im, _scc));

      _scc_repr[c] = SccRepr(i, st, sg);

      found_repr[c] = 1;
    }
  }
}

std::ostream &operator<<(std::ostream &os,
                         PartialPermInverseSemigroup const &ppisg)
{
  os << "GENERATORS: " << DUMP(ppisg.generators());
  return os;
}

} // namespace mpsym
