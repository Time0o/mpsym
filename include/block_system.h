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
  typedef std::vector<std::vector<unsigned>>::const_iterator const_iterator;

  BlockSystem() : _n(0u) {};
  BlockSystem(unsigned n, std::vector<std::vector<unsigned>> const &blocks);
  BlockSystem(std::vector<unsigned> const &classes);

  std::vector<unsigned> classes() const { return _classes; }
  std::vector<std::vector<unsigned>> blocks() const { return _blocks; }

  unsigned degree() const { return _n; }
  unsigned size() const { return static_cast<unsigned>(_blocks.size()); }
  bool trivial() const;

  std::vector<unsigned> const& operator[](unsigned const i) const;
  const_iterator begin() const;
  const_iterator end() const;

  PermGroup block_permuter(PermSet const &generators) const;

  static bool is_block(PermSet const &generators,
                       std::vector<unsigned> const &block);

  static std::vector<Perm> block_stabilizers(PermSet const &generators,
                                             std::vector<unsigned> const &block);

  static BlockSystem from_block(PermSet const &generators,
                                std::vector<unsigned> const &block);

  static BlockSystem minimal(PermSet const &generators,
                             std::vector<unsigned> const &initial_class);

  static std::vector<BlockSystem> non_trivial(PermGroup const &pg,
                                              bool assume_transitivity = false);

private:
  static std::vector<BlockSystem> non_trivial_transitive(PermGroup const &pg);
  static std::vector<BlockSystem> non_trivial_non_transitive(PermGroup const &pg);

  unsigned _n;
  std::vector<unsigned> _classes;
  std::vector<std::vector<unsigned>> _blocks;
};

std::ostream& operator<<(std::ostream& stream, BlockSystem const &bs);

} // namespace cgtl

#endif // _GUARD_BLOCK_SYSTEM_H
