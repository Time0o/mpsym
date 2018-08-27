#include <cassert>
#include <algorithm>
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

  bool transitive;

  if (assume_transitivity) {
    transitive = true;
    Dbg(Dbg::DBG) << "Assuming transitivity";
  } else {
    transitive = pg.transitive();
    Dbg(Dbg::DBG) << "Group " << (transitive ? " is " : " is not ") << " transitive";
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

  Dbg(Dbg::DBG) << "Finding all non-trivial block systems for group with generators:";
  Dbg(Dbg::DBG) << sgs;

  std::vector<Perm> stab1 = pg.bsgs().stabilizers(1);
  Dbg(Dbg::TRACE) << "G_1 has generators: " << stab1;

  std::vector<std::vector<unsigned>> stab1_orbits = SchreierSims::orbits(stab1);
  Dbg(Dbg::TRACE) << "Orbit decomposition of G_1 is: " << stab1_orbits;

  std::vector<BlockSystem> res;

  for (auto const &orbit : stab1_orbits) {
    unsigned repr = orbit[0];
    if (repr == 1u)
      continue;

    auto bs = BlockSystem::minimal(sgs, {1, repr});

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
  throw std::logic_error("not implemented");
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
