#ifndef _GUARD_BLOCK_SYSTEM_H
#define _GUARD_BLOCK_SYSTEM_H

#include <ostream>
#include <vector>

#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"

namespace cgtl
{

class BlockSystem
{
friend std::ostream& operator<<(std::ostream& stream, BlockSystem const &bs);

public:
  using Block = std::vector<unsigned>;

  typedef std::vector<Block>::const_iterator const_iterator;

  BlockSystem() : _n(0u) {};
  BlockSystem(unsigned n, std::vector<Block> const &blocks);
  BlockSystem(std::vector<unsigned> const &classes);

  std::vector<unsigned> classes() const { return _classes; }
  std::vector<Block> blocks() const { return _blocks; }

  unsigned degree() const { return _n; }
  unsigned size() const { return static_cast<unsigned>(_blocks.size()); }
  bool trivial() const;

  Block const& operator[](unsigned const i) const;
  const_iterator begin() const;
  const_iterator end() const;

  PermGroup block_permuter(PermSet const &generators) const;

  static bool is_block(PermSet const &generators, Block const &block);

  static PermSet block_stabilizers(PermSet const &generators,
                                   Block const &block);

  static BlockSystem from_block(PermSet const &generators,
                                Block const &block);

  static BlockSystem minimal(PermSet const &generators,
                             std::vector<unsigned> const &initial_class);

  static std::vector<BlockSystem> non_trivial(PermGroup const &pg,
                                              bool assume_transitivity = false);

private:
  static unsigned minimal_find_rep(unsigned k,
                              std::vector<unsigned> *classpath);

  static bool minimal_merge_classes(unsigned k1,
                                    unsigned k2,
                                    std::vector<unsigned> *classpath,
                                    std::vector<unsigned> *cardinalities,
                                    std::vector<unsigned> *queue);

  static std::vector<BlockSystem> non_trivial_transitive(PermGroup const &pg);

  static std::vector<BlockSystem> non_trivial_non_transitive(PermGroup const &pg);

  static Block non_trivial_shift_block(Block const &block, unsigned shift);

  static std::vector<Block> non_trivial_find_representatives(
    PermSet const &generators,
    std::vector<std::vector<BlockSystem>> const &partial_blocksystems,
    std::vector<unsigned> const &domain_offsets);

  static std::vector<BlockSystem> non_trivial_from_representatives(
    PermSet const &generators,
    std::vector<Block> const &representatives);

  unsigned _n;
  std::vector<unsigned> _classes;
  std::vector<std::vector<unsigned>> _blocks;
};

std::ostream& operator<<(std::ostream& stream, BlockSystem const &bs);

} // namespace cgtl

#endif // _GUARD_BLOCK_SYSTEM_H
