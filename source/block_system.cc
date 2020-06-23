#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

#include "block_system.h"
#include "dbg.h"
#include "dump.h"
#include "orbits.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"

namespace mpsym
{

bool BlockSystem::trivial() const
{ return size() == 1u || (*this)[0].size() == 1u; }

BlockSystem::Block const& BlockSystem::operator[](unsigned i) const
{ return _blocks[i]; }

BlockSystem::const_iterator BlockSystem::begin() const
{ return _blocks.begin(); }

BlockSystem::const_iterator BlockSystem::end() const
{ return _blocks.end(); }

unsigned BlockSystem::block_index(unsigned x) const
{ return _block_indices[x - 1u]; }

PermSet BlockSystem::block_permuter(PermSet const &generators_) const
{
  PermSet generators(generators_);

  std::vector<unsigned> perm(size());
  for (auto i = 0u; i < generators.size(); ++i) {
    auto gen(generators[i]);

    for (unsigned j = 0u; j < size(); ++j)
      perm[j] = block_index(gen[(*this)[j][0]]) + 1u;

    generators[i] = Perm(perm);
  }

  return generators;
}

PermSet BlockSystem::block_stabilizers(PermSet const &generators,
                                       Block const &block)
{
  // initialize block stabilizer generating set as generators of subgroup
  // stabilizing a block element (we arbitrarily choose the first one)
  PermGroup pg(generators.degree(), generators);

  pg.bsgs().base_change({block[0]});

  auto stabilizer_generators(pg.bsgs().stabilizers(0));
  auto stabilizer_orbit(Orbit::generate(block[0], stabilizer_generators));

  // extend block stabilizer generating set
  std::unordered_set<unsigned> block_elements(block.begin(), block.end());

  for (unsigned beta : pg.bsgs().orbit(0)) {
    if (block_elements.find(beta) == block_elements.end())
      continue;

    if (stabilizer_orbit.contains(beta))
      continue;

    Perm transv(pg.bsgs().transversal(0, beta));

    stabilizer_generators.insert(transv);
    stabilizer_orbit.update(stabilizer_generators, {transv});
  }

  return stabilizer_generators;
}

BlockSystem BlockSystem::minimal(PermSet const &generators,
                                 std::vector<unsigned> const &initial_block)
{
  assert(initial_block.size() >= 2u);

  std::vector<unsigned> classpath(generators.degree());
  std::vector<unsigned> cardinalities(generators.degree());
  std::vector<unsigned> queue;

  DBG(DEBUG) << "Finding minimal block system for:";
  DBG(DEBUG) << generators;

  DBG(TRACE) << "Initial block: " << initial_block;

  for (auto i = 0u; i < generators.degree(); ++i) {
    classpath[i] = i;
    cardinalities[i] = 1u;
  }

  for (auto i = 0u; i < initial_block.size() - 1u; ++i) {
    unsigned tmp = initial_block[i + 1u];

    classpath[tmp] = initial_block[0];
	queue.push_back(tmp);
  }

  cardinalities[initial_block[0]] = static_cast<unsigned>(initial_block.size());

  DBG(TRACE) << "Initial classpath: " << classpath;
  DBG(TRACE) << "Initial cardinalities: " << cardinalities;
  DBG(TRACE) << "Initial queue: " << queue;

  unsigned i = 0u;
  unsigned l = initial_block.size() - 2u;

  while (i <= l) {
    unsigned gamma = queue[i++];
    DBG(TRACE) << "Gamma: " << gamma;

    for (Perm const &gen : generators) {
      DBG(TRACE) << "Gen: " << gen;

      unsigned c1 = gen[gamma + 1u] - 1u;
      unsigned c2 = gen[minimal_find_rep(gamma, classpath) + 1u] - 1u;

      DBG(TRACE) << "Considering classes " << c1 << " and " << c2;

      if (minimal_merge_classes(c1, c2, classpath, cardinalities, queue))
        ++l;
    }
  }

  for (auto i = 0u; i < generators.degree(); ++i)
    minimal_find_rep(i, classpath);

  minimal_compress_classpath(classpath);

  DBG(TRACE) << "Final classpath is: " << classpath;

  BlockSystem res(classpath);

  DBG(DEBUG) << "=> Resulting minimal block system:";
  DBG(DEBUG) << res;

  return res;
}

std::vector<BlockSystem> BlockSystem::non_trivial(PermGroup const &pg,
                                                  bool assume_transitivity)
{
  assert((!assume_transitivity || pg.is_transitive()) &&
    "transitivity assumption correct");

  DBG(DEBUG) << "Finding all non-trivial block systems for:";
  DBG(DEBUG) << pg;

  bool transitive;

  if (assume_transitivity) {
    transitive = true;
    DBG(TRACE) << "Assuming transitivity";
  } else {
    transitive = pg.is_transitive();
    DBG(TRACE) << "Group " << (transitive ? "is" : "is not") << " transitive";
  }

  auto res(transitive ? non_trivial_transitive(pg)
                      : non_trivial_non_transitive(pg));

  DBG(DEBUG) << "=> Resulting non-trivial block systems:";
#ifndef NDEBUG
  for (BlockSystem const &bs : res)
    DBG(DEBUG) << bs;
#endif

  return res;
}

BlockSystem::BlockSystem(BlockIndices const &block_indices)
: _degree(block_indices.size()),
  _block_indices(block_indices)
{
  for (auto i = 1u; i <= degree(); ++i) {
    unsigned j = block_index(i);

    if (j + 1u > size()) {
      for (auto k = j + 1u - size(); k > 0u; --k)
        _blocks.emplace_back();
    }

    _blocks[j].push_back(i);
  }

  assert_blocks();
  assert_block_indices();
}

void BlockSystem::assert_blocks() const
{
#ifndef NDEBUG
  assert(size() > 0u && "number of blocks is positive");

  for (auto const &block : *this) {
    for (unsigned x : block)
      assert(x > 0u && "blocks have valid elements");
  }

  for (auto const &block : *this)
    assert(std::is_sorted(block.begin(), block.end()) && "blocks are sorted");

  for (auto i = 1u; i < size(); ++i)
    assert((*this)[i].size() == (*this)[0].size() && "blocks have same size");

  std::unordered_set<unsigned> block_union;
  for (auto const &block : *this)
    block_union.insert(block.begin(), block.end());

  assert(block_union.size() == degree() && "blocks partition domain");
#endif
}

void BlockSystem::assert_block_indices() const
{
  assert(_block_indices.size() == degree());

  for (unsigned x = 1u; x < degree(); ++x) {
    unsigned i = block_index(x);
    assert(i < size());

    auto block((*this)[i]);
    assert(std::find(block.begin(), block.end(), x) != block.end());
  }
}

bool BlockSystem::is_block(PermSet const &generators, Block const &block)
{
  auto maps_to_other_block = [&](Perm const &perm, unsigned x) {
    return std::find(block.begin(), block.end(), perm[x]) == block.end();
  };

  for (Perm const &gen : generators) {
    bool should_map_to_other_block = maps_to_other_block(gen, block[0]);

    for (auto i = 1u; i < block.size(); ++i) {
      if (maps_to_other_block(gen, block[i]) != should_map_to_other_block)
        return false;
    }
  }

  return true;
}

BlockSystem BlockSystem::from_block(PermSet const &generators,
                                    Block const &block)
{
  assert(is_block(generators, block));

  std::vector<Block> blocks{block};

  std::vector<int> block_indices(generators.degree(), -1);
  for (unsigned x : block)
    block_indices[x - 1u] = 0;

  unsigned current_block_idx = 0u;
  unsigned processed = block.size();

  while (current_block_idx < blocks.size()) {
    auto current_block(blocks[current_block_idx]);

    for (Perm const &gen : generators) {
      if (block_indices[gen[current_block[0]] - 1u] != -1)
        continue;

      blocks.push_back(current_block.permuted(gen));

      for (unsigned x : current_block)
        block_indices[gen[x] - 1u] = blocks.size();

      if ((processed += block.size()) == generators.degree())
        return BlockSystem(blocks.begin(), blocks.end());
    }

    ++current_block_idx;
  }

  throw std::logic_error("unreachable");
}

unsigned BlockSystem::minimal_find_rep(unsigned k,
                                       std::vector<unsigned> &classpath)
{
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
}

bool BlockSystem::minimal_merge_classes(unsigned k1,
                                        unsigned k2,
                                        std::vector<unsigned> &classpath,
                                        std::vector<unsigned> &cardinalities,
                                        std::vector<unsigned> &queue)
{
  unsigned r1 = minimal_find_rep(k1, classpath);
  unsigned r2 = minimal_find_rep(k2, classpath);

  DBG(TRACE) << "Representatives are: "
             << k1 << " => " << r1 << ", " << k2 << " => " << r2;

  if (r1 == r2)
    return false;

  if (cardinalities[r1] < cardinalities[r2])
    std::swap(r1, r2);

  DBG(TRACE) << "Merging classes:";

  classpath[r2] = r1;
  DBG(TRACE) << "Updated classpath: " << classpath;

  cardinalities[r1] += cardinalities[r2];
  DBG(TRACE) << "Updated cardinalities: " << cardinalities;

  queue.push_back(r2);
  DBG(TRACE) << "Updated queue: " << queue;

  return true;
}

void BlockSystem::minimal_compress_classpath(std::vector<unsigned> &classpath)
{
  std::unordered_map<unsigned, unsigned> compression;

  unsigned i = 0u;
  for (unsigned j : classpath) {
    if (compression.find(j) == compression.end())
      compression[j] = i++;
  }

  for (unsigned &j : classpath)
    j = compression[j];
}

std::vector<BlockSystem> BlockSystem::non_trivial_transitive(
  PermGroup const &pg)
{
  // first base element
  unsigned first_base_elem = pg.bsgs().base_point(0);
  DBG(TRACE) << "First base element is: " << first_base_elem;

  // generators of stabilizer subgroup for first base element
  PermSet stab = pg.bsgs().stabilizers(1);

  if (stab.empty()) {
    DBG(TRACE) << "No generators stabilizing first base element";
    return {};
  }

  DBG(TRACE) << "Generators stabilizing first base element:";
  DBG(TRACE) << stab;

  std::vector<BlockSystem> res;

  // iterate over orbit partition of stabilizer subgroup
  for (auto const &orbit : OrbitPartition(stab.degree(), stab)) {
    if (orbit[0] == first_base_elem)
      continue;

    // find minimal blocksystem corresponding to orbit
    auto bs(BlockSystem::minimal(pg.generators(), {first_base_elem, orbit[0]}));

    if (!bs.trivial()) {
      DBG(TRACE) << "Found blocksystem:";
      DBG(TRACE) << bs;
      res.push_back(bs);
    }
  }

  return res;
}

std::vector<BlockSystem> BlockSystem::non_trivial_non_transitive(
  PermGroup const &pg)
{
  OrbitPartition orbits(pg.degree(), pg.generators());

  DBG(TRACE) << "Group has " << orbits.num_partitions() << " distinct orbits:";
#ifndef NDEBUG
  for (auto const &orbit : orbits)
    DBG(TRACE) << orbit;
#endif

  std::vector<std::vector<BlockSystem>> partial_blocksystems(orbits.num_partitions());
  std::vector<unsigned> domain_offsets(orbits.num_partitions());

  for (auto i = 0u; i < orbits.num_partitions(); ++i) {
    // calculate all non trivial block systems for orbit restricted group
    PermSet restricted_gens;

    auto orbit_extremes =
      std::minmax_element(orbits[i].begin(), orbits[i].end());

    unsigned orbit_low = *std::get<0>(orbit_extremes);
    unsigned orbit_high = *std::get<1>(orbit_extremes);

    domain_offsets[i] = orbit_low - 1u;

    for (Perm const &gen : pg.generators()) {
      Perm perm(gen.restricted(orbits[i].begin(), orbits[i].end()));

      if (!perm.id())
        restricted_gens.insert(perm.normalized(orbit_low, orbit_high));
    }

    DBG(TRACE) << "Group generators restricted to " << orbits[i] << ":";
    DBG(TRACE) << restricted_gens;

    PermGroup pg_restricted(orbit_high - orbit_low + 1u, restricted_gens);
    partial_blocksystems[i] = non_trivial(pg_restricted, true);

    // append trivial blocksystem {{x} | x in orbit}
    std::vector<unsigned> trivial_classes(orbits[i].size());
    std::iota(trivial_classes.begin(), trivial_classes.end(), 0u);

    BlockSystem bs(trivial_classes);
    partial_blocksystems[i].push_back(bs);
  }

  DBG(TRACE) << "Relevant block systems for all group restrictions:";
#ifndef NDEBUG
  for (auto const &bs : partial_blocksystems)
    DBG(TRACE) << bs;
#endif

  auto blocksystem_representatives(
    non_trivial_find_representatives(pg.generators(),
                                     partial_blocksystems,
                                     domain_offsets));

  return non_trivial_from_representatives(pg.generators(),
                                          blocksystem_representatives);
}

std::vector<BlockSystem::Block> BlockSystem::non_trivial_find_representatives(
  PermSet const &generators,
  std::vector<std::vector<BlockSystem>> const &partial_blocksystems,
  std::vector<unsigned> const &domain_offsets)
{
  DBG(TRACE) << "Finding block system representatives";

  std::vector<Block> res;

  std::function<void(std::vector<BlockSystem const *> const &, unsigned, bool)>
  recurse = [&](std::vector<BlockSystem const *> const &current_blocksystems,
                unsigned i,
                bool one_trivial)
  {
    if (i == partial_blocksystems.size()) {
      DBG(TRACE) << "Considering block system combination:";
#ifndef NDEBUG
      for (auto j = 0u; j < current_blocksystems.size(); ++j) {
        DBG(TRACE) << *current_blocksystems[j]
                   << " (shifted by " << domain_offsets[j] << ")";
      }
#endif

      Block current_block_unshifted((*current_blocksystems[0])[0]);
      Block current_block(current_block_unshifted.shifted(domain_offsets[0]));

      for (auto j = 1u; j < current_blocksystems.size(); ++j) {
        bool next_block = false;

        for (auto const &block : *current_blocksystems[j]) {
          Block shifted_block(block.shifted(domain_offsets[j]));
          Block extended_block(current_block.unified(shifted_block));

          if (is_block(generators, extended_block)) {
            current_block = extended_block;
            next_block = true;
            break;
          }
        }

        if (!next_block)
          return;
      }

      res.push_back(current_block);
      DBG(TRACE) << "Found representative block: " << current_block;

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

  return res;
}

std::vector<BlockSystem> BlockSystem::non_trivial_from_representatives(
  PermSet const &generators,
  std::vector<Block> const &representatives)
{
  std::vector<BlockSystem> res;
  for (auto const &repr : representatives)
    res.push_back(from_block(generators, repr));

  return res;
}

std::ostream &operator<<(std::ostream &os, BlockSystem const &bs)
{
  os << DUMP_CUSTOM(bs._blocks, "{}", "{}");
  return os;
}

} // namespace mpsym
