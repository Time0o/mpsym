#ifndef _GUARD_BLOCK_SYSTEM_H
#define _GUARD_BLOCK_SYSTEM_H

#include <ostream>
#include <vector>

#include "perm.h"
#include "perm_group.h"

namespace cgtl
{

class BlockSystem
{
friend std::ostream& operator<<(std::ostream& stream, BlockSystem const &bs);

public:
  typedef std::vector<std::vector<unsigned>>::const_iterator const_iterator;

  BlockSystem(std::vector<unsigned> const &classes);

  unsigned degree() const { return _n; }
  unsigned size() const { return static_cast<unsigned>(_blocks.size()); }
  bool trivial() const { return _blocks.size() == 1u || _blocks[0].size() == 1u; }

  std::vector<unsigned> const& operator[](unsigned const i) const;
  const_iterator begin() const;
  const_iterator end() const;

  static bool is_block(std::vector<Perm> const &generators,
                       std::vector<unsigned> const &block);

  static BlockSystem minimal(std::vector<Perm> const &generators,
                             std::vector<unsigned> const &initial_class);

  static std::vector<BlockSystem> non_trivial(
    PermGroup const &pg, bool assume_transitivity = false);

private:
  static std::vector<BlockSystem> non_trivial_transitive(PermGroup const &pg);
  static std::vector<BlockSystem> non_trivial_non_transitive(PermGroup const &pg);

  unsigned _n;
  std::vector<std::vector<unsigned>> _blocks;
};

std::ostream& operator<<(std::ostream& stream, BlockSystem const &bs);

} // namespace cgtl

#endif // _GUARD_BLOCK_SYSTEM_H
