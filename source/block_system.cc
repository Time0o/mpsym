#include <algorithm>
#include <cassert>
#include <functional>
#include <utility>
#include <vector>

#include "block_system.h"
#include "dbg.h"
#include "orbits.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"

namespace cgtl
{

BlockSystem::BlockSystem(unsigned n, std::vector<Block> const &blocks)
: _n(n),
  _classes(n),
  _blocks(blocks)
{
#ifndef NDEBUG
  for (auto const &block : _blocks)
    assert(is_sorted(block.begin(), block.end()));
#endif

  for (auto i = 0u; i < blocks.size(); ++i) {
    for (unsigned x : blocks[i])
      _classes[x - 1u] = i + 1u;
  }
}

BlockSystem::BlockSystem(std::vector<unsigned> const &classes)
: _n(classes.size()),
  _classes(classes)
{
  std::vector<int> block_indices(_n, -1);

  for (auto i = 0u; i < _n; ++i) {
    unsigned c = _classes[i];

    if (block_indices[c - 1u] == -1) {
      _blocks.push_back({i + 1u});
      block_indices[c - 1u] = static_cast<int>(_blocks.size()) - 1;
    } else {
      _blocks[block_indices[c - 1u]].push_back(i + 1u);
    }
  }

#ifndef NDEBUG
  for (auto const &block : _blocks) {
    assert(block.size() == _blocks[0].size() &&
      "blocks in block system have same size");
  }
#endif
}

bool BlockSystem::trivial() const
{ return _blocks.size() == 1u || _blocks[0].size() == 1u; }

BlockSystem::Block const& BlockSystem::operator[](unsigned const i) const
{ return _blocks[i]; }

BlockSystem::const_iterator BlockSystem::begin() const
{ return _blocks.begin(); }

BlockSystem::const_iterator BlockSystem::end() const
{ return _blocks.end(); }

PermGroup BlockSystem::block_permuter(PermSet const &generators) const
{
  unsigned d = static_cast<unsigned>(_blocks.size());

  std::vector<unsigned> block_pointers(generators.degree() + 1u);

  for (unsigned i = 0u; i < d; ++i) {
    for (unsigned x : _blocks[i])
      block_pointers[x] = i + 1u;
  }

  std::vector<unsigned> tmp_perm(d);
  std::vector<Perm> permuter_generators(generators.size());

  for (auto i = 0u; i < generators.size(); ++i) {
    auto gen = generators[i];

    for (unsigned j = 0u; j < d; ++j) {
      unsigned x = _blocks[j][0];
      unsigned y = gen[x];

      tmp_perm[j] = block_pointers[y];
    }

    permuter_generators[i] = Perm(tmp_perm);
  }

  return PermGroup(
    d, PermSet(permuter_generators.begin(), permuter_generators.end()));
}

bool BlockSystem::is_block(PermSet const &generators, Block const &block)
{
  for (Perm const &gen : generators) {
    bool maps_to_other_block =
      std::find(block.begin(), block.end(), gen[block[0]]) == block.end();

    for (auto i = 0u; i < block.size(); ++i) {
      unsigned x = block[i];

      if (std::find(block.begin(), block.end(), gen[x]) == block.end()) {
        if (!maps_to_other_block)
          return false;
      } else {
        if (maps_to_other_block)
          return false;
      }
    }
  }

  return true;
}

PermSet BlockSystem::block_stabilizers(PermSet const &generators,
                                       Block const &block)
{
  throw std::logic_error(
    "TODO: block stabilizer group need to be found via backtracking");

  PermSet res;

  for (Perm const &gen : generators) {
    bool id = true;
    bool stabilizes = true;

    for (unsigned x : block) {
      unsigned y = gen[x];

      if (id && y != x)
        id = false;

      if (std::find(block.begin(), block.end(), y) == block.end()) {
        stabilizes = false;
        break;
      }
    }

    if (!id && stabilizes)
      res.insert(gen.restricted(block));
  }

  return res;
}

BlockSystem BlockSystem::from_block(PermSet const &generators,
                                    Block const &block)
{
  assert(is_block(generators, block));

  std::vector<unsigned> classes(generators.degree() + 1u, 0u);

  for (unsigned x : block)
    classes[x] = 1u;

  std::vector<Block> blocks {block};

  unsigned block_idx = 0u;
  unsigned n_processed = block.size();

  while (block_idx < blocks.size()) {
    auto current_block(blocks[block_idx]);

    for (Perm const &gen : generators) {
      unsigned x = blocks[block_idx][0];
      unsigned y = gen[x];

      if (classes[y] != 0u)
        continue;

      Block next_block(block.size());
      for (auto i = 0u; i < current_block.size(); ++i) {
        unsigned tmp = gen[current_block[i]];
        next_block[i] = tmp;
        classes[tmp] = blocks.size() + 1u;
      }
      blocks.push_back(next_block);

      if ((n_processed += block.size()) == generators.degree())
        return BlockSystem(generators.degree(), blocks);
    }

    ++block_idx;
  }

  throw std::logic_error("unreachable");
}

BlockSystem BlockSystem::minimal(PermSet const &generators,
                                 std::vector<unsigned> const &initial_class)
{
  assert(initial_class.size() >= 2u);

  std::vector<unsigned> classpath(generators.degree());
  std::vector<unsigned> cardinalities(generators.degree());
  std::vector<unsigned> queue;

  Dbg(Dbg::DBG) << "Finding minimal block system for:";
  Dbg(Dbg::DBG) << generators;

  Dbg(Dbg::DBG) << "Initial class is: " << initial_class;

  for (auto i = 1u; i <= generators.degree(); ++i) {
    classpath[i - 1u] = i;
    cardinalities[i - 1u] = 1u;
  }

  for (auto i = 0u; i < initial_class.size() - 1u; ++i) {
    unsigned tmp = initial_class[i + 1u];

    classpath[tmp - 1u] = initial_class[0];
	queue.push_back(tmp);
  }

  cardinalities[initial_class[0] - 1u] =
    static_cast<unsigned>(initial_class.size());

  Dbg(Dbg::TRACE) <<  "Initial classpath: " << classpath;
  Dbg(Dbg::TRACE) <<  "Initial cardinalities: " << cardinalities;
  Dbg(Dbg::TRACE) <<  "Initial queue: " << queue;

  unsigned i = 0u;
  unsigned l = initial_class.size() - 2u;

  while (i <= l) {
    unsigned gamma = queue[i++];
    Dbg(Dbg::TRACE) << "== Gamma: " << gamma;

    for (Perm const &gen : generators) {
      Dbg(Dbg::TRACE) << "= gen: " << gen;

      unsigned c1 = gen[gamma];
      unsigned c2 = gen[minimal_find_rep(gamma, &classpath)];

      Dbg(Dbg::TRACE) << "Considering classes " << c1 << " and " << c2;

      if (minimal_merge_classes(c1, c2, &classpath, &cardinalities, &queue))
        ++l;
    }
  }

  for (auto i = 1u; i <= generators.degree(); ++i)
    minimal_find_rep(i, &classpath);

  Dbg(Dbg::TRACE) << "Final classpath is: " << classpath;

  BlockSystem res(classpath);

  Dbg(Dbg::TRACE) << "==> Resulting minimal block system: " << res;

  return res;
}

std::vector<BlockSystem> BlockSystem::non_trivial(PermGroup const &pg,
                                                  bool assume_transitivity)
{
  assert((!assume_transitivity || pg.is_transitive()) &&
    "transitivity assumption correct");

  Dbg(Dbg::DBG) << "Finding all non-trivial block systems for:";
  Dbg(Dbg::DBG) << pg;

  bool transitive;

  if (assume_transitivity) {
    transitive = true;
    Dbg(Dbg::DBG) << "Assuming transitivity";
  } else {
    transitive = pg.is_transitive();
    Dbg(Dbg::DBG) << "Group " << (transitive ? "is" : "is not") << " transitive";
  }

  if (transitive) {
    return non_trivial_transitive(pg);
  } else {
    return non_trivial_non_transitive(pg);
  }
}

unsigned BlockSystem::minimal_find_rep(unsigned k,
                                       std::vector<unsigned> *classpath)
{
  // find class
  unsigned res = k;
  unsigned next = (*classpath)[res - 1u];

  while (next != res) {
    res = next;
    next = (*classpath)[res - 1u];
  }

  // compress path
  unsigned current = k;
  next = (*classpath)[k - 1u];

  while (next != current) {
    (*classpath)[current - 1u] = res;

    current = next;
    next = (*classpath)[current - 1u];
  }

  return res;
}

bool BlockSystem::minimal_merge_classes(unsigned k1,
                                        unsigned k2,
                                        std::vector<unsigned> *classpath,
                                        std::vector<unsigned> *cardinalities,
                                        std::vector<unsigned> *queue)
{
  unsigned r1 = minimal_find_rep(k1, classpath);
  unsigned r2 = minimal_find_rep(k2, classpath);

  Dbg(Dbg::TRACE) << "Representatives are: "
                  << k1 << " => " << r1 << ", " << k2 << " => " << r2;

  if (r1 == r2)
    return false;

  if ((*cardinalities)[r1 - 1u] < (*cardinalities)[r2 - 1u])
    std::swap(r1, r2);

  Dbg(Dbg::TRACE) << "=> Merging classes:";

  (*classpath)[r2 - 1u] = r1;
  Dbg(Dbg::TRACE) << "Updated classpath: " << *classpath;

  (*cardinalities)[r1 - 1u] += (*cardinalities)[r2 - 1u];
  Dbg(Dbg::TRACE) << "Updated cardinalities: " << *cardinalities;

  queue->push_back(r2);
  Dbg(Dbg::TRACE) << "Updated queue: " << *queue;

  return true;
}

std::vector<BlockSystem> BlockSystem::non_trivial_transitive(
  PermGroup const &pg)
{
  // TODO: does the first base element HAVE to be one?
  unsigned first_base_elem = pg.bsgs().base_point(0);
  Dbg(Dbg::TRACE) << "First base element is: " << first_base_elem;

  PermSet stab = pg.bsgs().stabilizers(1);

  if (stab.empty()) {
    Dbg(Dbg::TRACE)
      << "No generators stabilizing first base element, no block systems possible";

    return std::vector<BlockSystem>();
  }

  Dbg(Dbg::TRACE) << "Generators stabilizing first base element: " << stab;

  std::vector<std::vector<unsigned>> stab_orbits = orbit_partition_expanded(stab);
  Dbg(Dbg::TRACE) << "Orbit decomposition of associated group is: " << stab_orbits;

  std::vector<BlockSystem> res;

  for (auto const &orbit : stab_orbits) {
    unsigned repr = orbit[0];
    if (repr == first_base_elem)
      continue;

    auto bs = BlockSystem::minimal(pg.bsgs().strong_generators(),
                                   {first_base_elem, repr});

    if (!bs.trivial()) {
      Dbg(Dbg::TRACE) << "Found blocksystem: " << bs;
      res.push_back(bs);
    }
  }

  Dbg(Dbg::TRACE) << "==> Resulting non-trivial block systems:";
  Dbg(Dbg::TRACE) << res;

  return res;
}

std::vector<BlockSystem> BlockSystem::non_trivial_non_transitive(
  PermGroup const &pg)
{
  auto orbits(orbit_partition_expanded(pg.bsgs().strong_generators()));
  auto gens(pg.bsgs().strong_generators());

  Dbg(Dbg::TRACE) << "Group has " << orbits.size() << " distinct orbits:";
#ifndef NDEBUG
  for (auto const &orbit : orbits)
    Dbg(Dbg::TRACE) << orbit;
#endif

  std::vector<std::vector<BlockSystem>> partial_blocksystems(orbits.size());
  std::vector<unsigned> domain_offsets(orbits.size());

  for (auto i = 0u; i < orbits.size(); ++i) {
    // calculate all non trivial block systems for orbit restricted group
    std::vector<Perm> restricted_gens;

    auto orbit_extremes =
      std::minmax_element(orbits[i].begin(), orbits[i].end());

    unsigned orbit_low = *std::get<0>(orbit_extremes);
    unsigned orbit_high = *std::get<1>(orbit_extremes);

    domain_offsets[i] = orbit_low - 1u;

    for (auto j = 0u; j < gens.size(); ++j) {
      Perm tmp(gens[j].restricted(orbits[i]));

      if (!tmp.id())
        restricted_gens.push_back(tmp.normalized(orbit_low, orbit_high));
    }

    Dbg(Dbg::TRACE) << "Group generators restricted to " << orbits[i] << ":";
	Dbg(Dbg::TRACE) << restricted_gens;

    PermGroup pg_restricted(
      orbit_high - orbit_low + 1u,
      PermSet(restricted_gens.begin(), restricted_gens.end()));

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

  Dbg(Dbg::TRACE) << "=> Relevant block systems for all group restrictions:";
#ifndef NDEBUG
  for (auto const &bs : partial_blocksystems)
    Dbg(Dbg::TRACE) << bs;
#endif

  auto representatives(
    non_trivial_find_representatives(gens, partial_blocksystems, domain_offsets));

  return non_trivial_from_representatives(gens, representatives);
}

BlockSystem::Block BlockSystem::non_trivial_shift_block(Block const &block,
                                                        unsigned shift)
{
  Block res(block.size());

  for (auto i = 0u; i < block.size(); ++i)
    res[i] = block[i] + shift;

  return res;
}

std::vector<BlockSystem::Block> BlockSystem::non_trivial_find_representatives(
  PermSet const &generators,
  std::vector<std::vector<BlockSystem>> const &partial_blocksystems,
  std::vector<unsigned> const &domain_offsets)
{
  Dbg(Dbg::TRACE) << "== Finding block system representatives";

  std::vector<Block> res;

  std::function<void(std::vector<BlockSystem const *> const &, unsigned, bool)>
  recurse = [&](std::vector<BlockSystem const *> const &current_blocksystems,
                unsigned i,
                bool one_trivial)
  {
    if (i == partial_blocksystems.size()) {
      Dbg(Dbg::TRACE) << "= Considering block system combination:";
#ifndef NDEBUG
      for (auto j = 0u; j < current_blocksystems.size(); ++j) {
        Dbg(Dbg::TRACE) << *current_blocksystems[j]
                        << " (shifted by " << domain_offsets[j] << ")";
      }
#endif

      Block current_block(non_trivial_shift_block((*current_blocksystems[0])[0],
                          domain_offsets[0]));

      for (auto j = 1u; j < current_blocksystems.size(); ++j) {
        bool next_block = false;

        BlockSystem const *bs = current_blocksystems[j];
        for (auto const &block : *bs) {
          Block shifted_block(non_trivial_shift_block(block, domain_offsets[j]));

          std::vector<unsigned> extended_block(
            current_block.size() + shifted_block.size());

          auto it = std::set_union(current_block.begin(), current_block.end(),
                                   shifted_block.begin(), shifted_block.end(),
                                   extended_block.begin());

          extended_block.resize(it - extended_block.begin());

          if (is_block(generators, extended_block)) {
            Dbg(Dbg::TRACE) << extended_block << " is a block";
            current_block = extended_block;
            next_block = true;
            break;
          }

          Dbg(Dbg::TRACE) << extended_block << " is not a block";
        }

        if (!next_block)
          return;
      }

      res.push_back(current_block);
      Dbg(Dbg::TRACE) << "=> Found representative block: " << current_block;

      return;
    }

    for (BlockSystem const &blocksystem : partial_blocksystems[i]) {
      if (blocksystem.trivial()) {
        if (one_trivial)
          return;

        one_trivial = true;
      }

      auto extended_blocksystems(current_blocksystems);
      extended_blocksystems.push_back(&blocksystem);

      recurse(extended_blocksystems, i + 1, one_trivial);
    }
  };

  recurse({}, 0u, false);

  Dbg(Dbg::TRACE) << "=> Representative blocks found:";
#ifndef NDEBUG
  for (auto const &block : res)
    Dbg(Dbg::TRACE) << block;
#endif

  return res;
}

std::vector<BlockSystem> BlockSystem::non_trivial_from_representatives(
  PermSet const &generators,
  std::vector<Block> const &representatives)
{
  Dbg(Dbg::TRACE) << "== Finding block systems from representatives";

  std::vector<BlockSystem> res(representatives.size());
  for (auto i = 0u; i < representatives.size(); ++i)
    res[i] = from_block(generators, representatives[i]);

  Dbg(Dbg::TRACE) << "==> Resulting block systems:";
#ifndef NDEBUG
  for (BlockSystem const &bs : res)
    Dbg(Dbg::TRACE) << bs;
#endif

  return res;
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
