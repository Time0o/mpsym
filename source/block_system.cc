#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "block_system.h"
#include "dbg.h"
#include "perm.h"
#include "perm_group.h"
#include "schreier_sims.h"

namespace cgtl
{

BlockSystem::BlockSystem(std::vector<unsigned> classes) : _n(classes.size())
{
  std::vector<int> block_indices(_n, -1);

  for (auto i = 0u; i < _n; ++i) {
    unsigned c = classes[i];

    if (block_indices[c - 1u] == -1) {
      _blocks.push_back({i + 1u});
      block_indices[c - 1u] = static_cast<int>(_blocks.size()) - 1;
    } else {
      _blocks[block_indices[c - 1u]].push_back(i + 1u);
    }
  }

  for (auto const &block : _blocks) {
    assert(block.size() == _blocks[0].size() &&
      "blocks in block system have same size");
  }
}

std::vector<unsigned> const& BlockSystem::operator[](unsigned const i) const
{ return _blocks[i]; }

BlockSystem::const_iterator BlockSystem::begin() const
{ return _blocks.begin(); }

BlockSystem::const_iterator BlockSystem::end() const
{ return _blocks.end(); }

BlockSystem BlockSystem::minimal(std::vector<Perm> const &generators,
                                 std::vector<unsigned> const &initial_class)
{
  assert(initial_class.size() >= 2u);

  unsigned degree = generators[0].degree();

  std::vector<unsigned> classpath(degree + 1u);
  std::vector<unsigned> cardinalities(degree + 1u);
  std::vector<unsigned> queue;

  auto rep = [&](unsigned k) {
    // find class
    unsigned res = k;
    unsigned next = classpath[res];

    while (next != res) {
      res = next;
      next = classpath[res];
    }

    // compress path
    unsigned current = k;
    next = classpath[k];
    while (next != current) {
      classpath[current] = res;
      current = next;
      next = classpath[current];
    }

    return res;
  };

  auto merge = [&](unsigned k1, unsigned k2) {
    unsigned r1 = rep(k1);
    unsigned r2 = rep(k2);

    Dbg(Dbg::TRACE) << "Representatives are: "
                    << k1 << " => " << r1 << ", " << k2 << " => " << r2;

    if (r1 == r2)
      return false;

    unsigned tmp1, tmp2;
    if (cardinalities[r1] >= cardinalities[r2]) {
      tmp1 = r1;
      tmp2 = r2;
    } else {
      tmp2 = r1;
      tmp1 = r2;
    }

    Dbg(Dbg::TRACE) << "=> Merging classes:";

    classpath[tmp2] = tmp1;
    Dbg(Dbg::TRACE) << "Updated classpath: "
                    << std::vector<unsigned>(classpath.begin() + 1,
                                             classpath.end());

    cardinalities[tmp1] += cardinalities[tmp2];
    Dbg(Dbg::TRACE) << "Updated cardinalities: "
                    << std::vector<unsigned>(cardinalities.begin() + 1,
                                             cardinalities.end());

	queue.push_back(tmp2);
    Dbg(Dbg::TRACE) << "Updated queue: " << queue;

    return true;
  };

  Dbg(Dbg::DBG) << "Finding minimal block system for:";
  Dbg(Dbg::DBG) << generators;

  Dbg(Dbg::DBG) << "Initial class is: " << initial_class;

  for (auto i = 1u; i <= degree; ++i) {
    classpath[i] = i;
    cardinalities[i] = 1u;
  }

  for (auto i = 0u; i < initial_class.size() - 1u; ++i) {
    unsigned tmp = initial_class[i + 1u];

    classpath[tmp] = initial_class[0];
	queue.push_back(tmp);
  }

  cardinalities[initial_class[0]] = static_cast<unsigned>(initial_class.size());

  Dbg(Dbg::TRACE) <<  "Initial classpath: "
                  << std::vector<unsigned>(classpath.begin() + 1,
                                           classpath.end());

  Dbg(Dbg::TRACE) <<  "Initial cardinalities: "
                  << std::vector<unsigned>(cardinalities.begin() + 1,
                                           cardinalities.end());

  Dbg(Dbg::TRACE) <<  "Initial queue: " << queue;

  unsigned i = 0u;
  unsigned l = initial_class.size() - 2u;

  while (i <= l) {
    unsigned gamma = queue[i++];
    Dbg(Dbg::TRACE) << "== Gamma: " << gamma;

    for (Perm const &gen : generators) {
      Dbg(Dbg::TRACE) << "= gen: " << gen;

      unsigned c1 = gen[gamma];
      unsigned c2 = gen[rep(gamma)];

      Dbg(Dbg::TRACE) << "Considering classes " << c1 << " and " << c2;

      if (merge(c1, c2))
        ++l;
    }
  }

  BlockSystem res(std::vector<unsigned>(classpath.begin() + 1,
                                        classpath.end()));

  Dbg(Dbg::TRACE) << "Resulting block system: " << res;

  return res;
}

std::vector<BlockSystem> BlockSystem::non_trivial(
  PermGroup const &pg, bool assume_transitivity)
{
  assert((!assume_transitivity || pg.transitive()) &&
    "transitivity assumption correct");

  Dbg(Dbg::DBG) << "Finding all non-trivial block systems for:";
  Dbg(Dbg::DBG) << pg;

  bool transitive;

  if (assume_transitivity) {
    transitive = true;
    Dbg(Dbg::DBG) << "Assuming transitivity";
  } else {
    transitive = pg.transitive();
    Dbg(Dbg::DBG) << "Group " << (transitive ? "is" : "is not") << " transitive";
  }

  if (transitive) {
    return non_trivial_transitive(pg);
  } else {
    return non_trivial_non_transitive(pg);
  }
}

std::vector<BlockSystem> BlockSystem::non_trivial_transitive(
  PermGroup const &pg)
{
  auto sgs = pg.bsgs().sgs();

  // TODO: does the first base element HAVE to be one?
  unsigned first_base_elem = pg.bsgs()[0].elem();
  Dbg(Dbg::TRACE) << "First base element is: " << first_base_elem;

  std::vector<Perm> stab = pg.bsgs().stabilizers(1);

  if (stab.empty()) {
    Dbg(Dbg::TRACE)
      << "No generators stabilizing first base element, no block systems possible";

    return std::vector<BlockSystem>();
  }

  Dbg(Dbg::TRACE) << "Generators stabilizing first base element: " << stab;

  std::vector<std::vector<unsigned>> stab_orbits = SchreierSims::orbits(stab);
  Dbg(Dbg::TRACE) << "Orbit decomposition of associated group is: " << stab_orbits;

  std::vector<BlockSystem> res;

  for (auto const &orbit : stab_orbits) {
    unsigned repr = orbit[0];
    if (repr == first_base_elem)
      continue;

    auto bs = BlockSystem::minimal(sgs, {first_base_elem, repr});

    if (!bs.trivial()) {
      Dbg(Dbg::TRACE) << "Found blocksystem: " << bs;
      res.push_back(bs);
    }
  }

  Dbg(Dbg::TRACE) << "Resulting block systems are:";
  Dbg(Dbg::TRACE) << res;

  return res;
}

std::vector<BlockSystem> BlockSystem::non_trivial_non_transitive(
  PermGroup const &pg)
{
  auto orbits(pg.orbits());
  auto gens(pg.bsgs().sgs());

  Dbg(Dbg::TRACE) << "Group has " << orbits.size() << " distinct orbits:";
#ifndef NDEBUG
  for (auto const &orbit : orbits)
    Dbg(Dbg::TRACE) << orbit;
#endif

  std::vector<std::vector<BlockSystem>> partial_blocksystems(orbits.size());

  for (auto i = 0u; i < orbits.size(); ++i) {
    // calculate all non trivial block systems for orbit restricted group
    std::vector<Perm> restricted_gens;

    auto orbit_extremes =
      std::minmax_element(orbits[i].begin(), orbits[i].end());

    unsigned orbit_low = *std::get<0>(orbit_extremes);
    unsigned orbit_high = *std::get<1>(orbit_extremes);

    for (auto j = 0u; j < gens.size(); ++j) {
      bool id;
      Perm tmp(gens[j].restricted(orbits[i], &id));

      if (!id)
        restricted_gens.push_back(tmp.shifted(orbit_low, orbit_high));
    }

    Dbg(Dbg::TRACE) << "Group generators restricted to " << orbits[i] << ":";
	Dbg(Dbg::TRACE) << restricted_gens;

    auto pg_restricted(PermGroup(orbit_high - orbit_low + 1u, restricted_gens));
    auto block_systems(non_trivial(pg_restricted, true));

    partial_blocksystems[i] = block_systems;

    Dbg(Dbg::TRACE) << "=> Resulting non-trivial block systems:";
#ifndef NDEBUG
    if (partial_blocksystems[i].size() == 0u) {
      Dbg(Dbg::TRACE) << "None";
    } else {
      for (auto const &bs : partial_blocksystems[i])
        Dbg(Dbg::TRACE) << bs;
   }
#endif

    // append trivial blocksystem {{x} | x in orbit}
    std::vector<unsigned> trivial_classes(orbits[i].size());
    for (auto j = 1u; j <= orbits[i].size(); ++j)
      trivial_classes[j - 1u] = j;

    partial_blocksystems[i].push_back(BlockSystem(trivial_classes));
  }

  Dbg(Dbg::TRACE) << "==> Relevant block systems for all group restrictions:";
#ifndef NDEBUG
  for (auto const &bs : partial_blocksystems)
    Dbg(Dbg::TRACE) << bs;
#endif

  return std::vector<BlockSystem>();
}

std::ostream& operator<<(std::ostream& stream, BlockSystem const &bs)
{
  stream << '{';
  for (auto i = 0u; i < bs._blocks.size(); ++i) {
    stream << '{' << bs._blocks[i][0];
    for (auto j = 1u; j < bs._blocks[i].size(); ++j)
      stream << ", " << bs._blocks[i][j];
    stream << '}';

    if (i != bs._blocks.size() - 1u)
      stream << ", ";
  }
  stream << '}';

  return stream;
}

} // namespace cgtl
