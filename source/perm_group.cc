#include <cassert>
#include <ctime>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "dbg.h"
#include "perm.h"

namespace cgtl
{

std::vector<unsigned> SchreierTree::orbit() const
{
  std::vector<unsigned> orbit {_root};

  for (auto const &node : _edges)
    orbit.push_back(node.first);

  return orbit;
}

Perm SchreierTree::transversal(unsigned origin) const
{
  Perm result(_degree);

  unsigned current = origin;
  while(current != _root) {
    Perm const &label = _labels.find(current)->second;
    result = label * result;
    current = _edges.find(current)->second;
  }

  return result;
}

std::vector<Perm> SchreierTree::transversals(
  std::vector<unsigned> const &origins) const
{
  std::vector<Perm> result;
  for (unsigned o : origins)
    result.push_back(transversal(o));

  return result;
}

PermGroup::PermGroup(unsigned degree, std::vector<Perm> const &generators)
  : _n(degree), _strong_generating_set(generators)
{
  schreier_sims(_base, _strong_generating_set, _schreier_trees);

  _order = 1u;
  for (auto const &st : _schreier_trees)
    _order *= st.orbit().size();
}

PermGroup PermGroup::symmetric(unsigned degree)
{
  assert(degree > 0u);

  if (degree == 1u)
    return PermGroup(1u, {Perm(1u)});

  std::vector<unsigned> gen;
  for (unsigned i = 1u; i <= degree; ++i)
    gen.push_back(i);

  return PermGroup(degree, {Perm(degree, {{1, 2}}), Perm(degree, {gen})});
}

PermGroup PermGroup::cyclic(unsigned degree)
{
  assert(degree > 0u);

  std::vector<unsigned> gen;
  for (unsigned i = 1u; i <= degree; ++i)
    gen.push_back(i);

  return PermGroup(degree, {Perm(degree, {gen})});
}

PermGroup PermGroup::alternating(unsigned degree)
{
  assert(degree > 2u);

  std::vector<Perm> gens;
  for (unsigned i = 3u; i <= degree; ++i)
    gens.push_back(Perm(degree, {{1, 2, i}}));

  return PermGroup(degree, gens);
}

std::vector<unsigned> PermGroup::orbit(unsigned alpha,
  std::vector<Perm> const &generators, SchreierTree &st)
{
  assert(generators.size() > 0u && "generator set not empty");
  assert(alpha <= generators[0].degree() && "alpha <= N");

  std::vector<unsigned> result {alpha};
  st.create_root(alpha);

  std::vector<unsigned> stack {alpha};
  std::set<unsigned> done {alpha};

  while (!stack.empty()) {
    unsigned beta = stack.back();
    stack.pop_back();

    for (auto i = 0u; i < generators.size(); ++i) {
      Perm const &gen = generators[i];
      unsigned beta_prime = gen[beta];

      if (done.find(beta_prime) == done.end()) {
        result.push_back(beta_prime);
        done.insert(beta_prime);
        stack.push_back(beta_prime);

        st.create_edge(beta_prime, beta, generators[i]);
      }
    }
  }

  return result;
}

std::pair<Perm, unsigned> PermGroup::strip(Perm const &perm,
  std::vector<unsigned> const &base, std::vector<Perm> const &generators,
  std::vector<SchreierTree> const &sts)
{
  assert(base.size() > 0u && "base not empty");
  assert(generators.size() > 0u && "generator set not empty");
  assert(sts.size() > 0u && "schreier tree set not empty");

  Perm result(perm);

  for (unsigned i = 0u; i < base.size(); ++i) {
    unsigned beta = result[base[i]];
    if (!sts[i].contains(beta))
      return std::make_pair(result, i + 1u);

    result *= ~sts[i].transversal(beta);
  }

  return std::make_pair(result, base.size() + 1u);
}

void PermGroup::schreier_sims(std::vector<unsigned> &base,
  std::vector<Perm> &generators, std::vector<SchreierTree> &schreier_trees)
{
  assert(generators.size() > 0u && "generator set not empty");

  Dbg(Dbg::DBG) << "Executing schreier sims algorithm";
  Dbg(Dbg::DBG) << "=== Input";
  Dbg(Dbg::DBG) << "B = " << base;
  Dbg(Dbg::DBG) << "S = " << generators;

  unsigned degree = generators[0].degree();

  // find initial base point if given base is empty
  if (base.size() == 0u) {
    bool init = false;
    for (unsigned beta = 1u; beta <= degree; ++beta) {
      for (auto const &gen : generators) {
        if (gen[beta] != beta) {
          base.push_back(beta);

          init = true;
          break;
        }
      }
      if (init)
        break;
    }
  }

  std::vector<std::vector<Perm>> strong_generators(base.size());
  std::vector<std::vector<unsigned>> fundamental_orbits(base.size());
  schreier_trees.resize(base.size(), SchreierTree(degree));

  // calculate initial strong generator sets
  for (unsigned i = 0u; i < base.size(); ++i) {
    for (Perm const &gen : generators) {
      bool stabilizes = true;
      for (unsigned k = 0u; k < i; ++k) {
        if (gen[base[k]] != base[k]) {
          stabilizes = false;
          break;
        }
      }
      if (stabilizes) {
        strong_generators[i].push_back(gen);
      }
    }

    fundamental_orbits[i] = orbit(
      base[i], strong_generators[i], schreier_trees[i]);
  }

  Dbg(Dbg::DBG) << "=== Initial values";
  Dbg(Dbg::DBG) << "B = " << base;
  for (unsigned i = 0u; i < base.size(); ++i) {
    Dbg(Dbg::DBG) << "S(" << (i + 1u) << ") = " << strong_generators[i];
    Dbg(Dbg::DBG) << "O(" << (i + 1u) << ") = " << fundamental_orbits[i];
    Dbg(Dbg::DBG) << "U(" << (i + 1u) << ") = "
                  << schreier_trees[i].transversals(fundamental_orbits[i]);
  }

  unsigned i = base.size();
  while (i >= 1u) {
top:
    Dbg(Dbg::DBG) << "=== Main loop (i = " << i << ")";
    auto const &st = schreier_trees[i - 1u];

    for (unsigned beta : fundamental_orbits[i - 1u]) {
      Dbg(Dbg::DBG) << "== Orbit element " << beta;

      Perm u_beta = st.transversal(beta);

      for (Perm const &x : strong_generators[i - 1u]) {
        unsigned beta_x = x[beta];

        Perm schreier_generator = u_beta * x * ~st.transversal(beta_x);
        Dbg(Dbg::DBG) << "Schreier Generator: " << schreier_generator;

        if (schreier_generator.id())
          continue;

        bool update_strong_generators = false;
        std::pair<Perm, unsigned> strip_result = strip(schreier_generator,
          base, generators, schreier_trees);

        Perm strip_perm = std::get<0>(strip_result);
        unsigned strip_level = std::get<1>(strip_result);
        Dbg(Dbg::DBG) << "Strips to: " << strip_perm << ", " << strip_level;

        if (strip_level <= base.size()) {
           update_strong_generators = true;

        } else if (strip_perm != Perm(degree)) {
          update_strong_generators = true;

          unsigned next_basepoint = 1u;
          while (strip_perm[next_basepoint] == next_basepoint)
            ++next_basepoint;

          base.push_back(next_basepoint);

          Dbg(Dbg::DBG) << "Adjoined new basepoint:";
          Dbg(Dbg::DBG) << "B = " << base;
        }

        if (update_strong_generators) {
          Dbg(Dbg::DBG) << "Updating strong generators:";

          for (unsigned j = i; j < strip_level; ++ j) {
            if (strong_generators.size() <= j) {
              strong_generators.push_back({});
              fundamental_orbits.push_back({});
              schreier_trees.push_back(SchreierTree(degree));
            }

            strong_generators[j].push_back(strip_perm);

            fundamental_orbits[j] =
              orbit(base[j], strong_generators[j], schreier_trees[j]);

            Dbg(Dbg::DBG) << "S(" << (j + 1u) << ") = " << strong_generators[j];
            Dbg(Dbg::DBG) << "O(" << (j + 1u) << ") = " << fundamental_orbits[j];
            Dbg(Dbg::DBG) << "U(" << (j + 1u) << ") = "
                          << schreier_trees[j].transversals(fundamental_orbits[j]);
          }

          i = strip_level;
          goto top;
        }
      }
    }
    i -= 1u;
  }

  generators.clear();
  for (auto const &gens : strong_generators)
    generators.insert(generators.end(), gens.begin(), gens.end());
}

bool PermGroup::is_element(Perm const &perm)
{
  Dbg(Dbg::DBG) << "Performing membership test for " << perm << " in:";
  Dbg(Dbg::DBG) << (*this);

  auto strip_result =
    strip(perm, _base, _strong_generating_set, _schreier_trees);

  Dbg(Dbg::TRACE) << "Strip returned " << std::get<0>(strip_result) << ", "
                  << std::get<1>(strip_result);

  bool ret = (std::get<1>(strip_result) == _base.size() + 1) &&
             (std::get<0>(strip_result).id());

  Dbg(Dbg::DBG) << (ret ? "=> Member" : "=> No Member");

  return ret;
}

Perm PermGroup::random_element()
{
  static std::default_random_engine gen(time(0));

  Perm result(_n);
  for (auto const &st : _schreier_trees) {
    std::vector<unsigned> orbit = st.orbit();
    std::uniform_int_distribution<> d(0u, orbit.size() - 1u);
    result *= st.transversal(orbit[d(gen)]);
  }

  return result;
}

PermGroup::const_iterator::const_iterator(
  std::vector<SchreierTree> schreier_trees) : _end(false)
{
  for (auto const &st : schreier_trees) {
    _state.push_back(0u);
    _transversals.push_back(st.transversals(st.orbit()));
    _current_factors.push_back(_transversals.back()[0]);
  }

  update_result();
}

PermGroup::const_iterator PermGroup::const_iterator::operator++()
{
  PermGroup::const_iterator pre(*this);
  next_state();
  return pre;
}

bool PermGroup::const_iterator::operator==(
  PermGroup::const_iterator const &rhs) const
{
  if (_end != rhs._end)
    return false;

  if (_end && rhs._end)
    return true;

  for (unsigned i = 0u; i < _state.size(); ++i) {
    if (_state[i] != rhs._state[i])
      return false;
  }
  return true;
}

void PermGroup::const_iterator::next_state()
{
  for (unsigned i = 0u; i < _state.size(); ++i) {
    _state[i]++;
    if (_state[i] == _transversals[i].size())
      _state[i] = 0u;

    _current_factors[i] = _transversals[i][_state[i]];

    if (i == _state.size() - 1u && _state[i] == 0u) {
      _end = true;
      break;
    }

    if (_state[i] != 0u)
      break;
  }

  update_result();
}

void PermGroup::const_iterator::update_result()
{
  _current_result = _current_factors[0];
  for (unsigned j = 1u; j < _current_factors.size(); ++j)
    _current_result *= _current_factors[j];
}

} // namespace cgtl
