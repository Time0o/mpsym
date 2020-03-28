#ifndef _GUARD_BLOCK_SYSTEM_H
#define _GUARD_BLOCK_SYSTEM_H

#include <initializer_list>
#include <ostream>
#include <vector>

#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"

namespace mpsym
{

class BlockSystem
{
  friend std::ostream &operator<<(std::ostream &os, BlockSystem const &bs);

public:
  struct Block : public std::vector<unsigned>
  {
    Block()
    : std::vector<unsigned>()
    {}

    Block(std::initializer_list<unsigned> block)
    : std::vector<unsigned>(block)
    {}

    Block permuted(Perm const &perm) const
    {
      Block res;
      for (unsigned x : *this)
        res.push_back(perm[x]);

      return res;
    }

    Block shifted(unsigned shift) const
    {
      Block res;
      for (unsigned x : *this)
        res.push_back(x + shift);

      return res;
    }

    Block unified(Block const &other) const
    {
      Block res;
      res.resize(size() + other.size());

      auto it = std::set_union(begin(), end(),
                               other.begin(), other.end(),
                               res.begin());

      res.resize(it - res.begin());

      return res;
    }
  };

  using BlockIndices = std::vector<unsigned>;

  using const_iterator = std::vector<Block>::const_iterator;

  unsigned degree() const { return _degree; }
  unsigned size() const { return static_cast<unsigned>(_blocks.size()); }
  bool trivial() const;

  Block const& operator[](unsigned i) const;
  const_iterator begin() const;
  const_iterator end() const;

  unsigned block_index(unsigned x) const;

  PermSet block_permuter(PermSet const &generators) const;

  static PermSet block_stabilizers(PermSet const &generators,
                                     Block const &block);

  static BlockSystem minimal(PermSet const &generators,
                             std::vector<unsigned> const &initial_block);

  static std::vector<BlockSystem> non_trivial(PermGroup const &pg,
                                              bool assume_transitivity = false);

private:
  template<typename IT>
  BlockSystem(IT first, IT last)
  : _blocks(first, last)
  {
    _degree = 0u;
    for (IT it = first; it != last; ++it) {
      for (unsigned x : *it)
        _degree = std::max(_degree, x);
    }

    _block_indices = std::vector<unsigned>(_degree);
    for (unsigned i = 0u; i < size(); ++i) {
      for (unsigned x : (*this)[i])
        _block_indices[x - 1u] = i;
    }

    assert_blocks();
    assert_block_indices();
  }

  BlockSystem(BlockIndices const &block_indices);

  void assert_blocks() const;
  void assert_block_indices() const;

  static bool is_block(PermSet const &generators, Block const &block);

  static BlockSystem from_block(PermSet const &generators, Block const &block);

  static unsigned minimal_find_rep(unsigned k,
                                   std::vector<unsigned> &classpath);

  static bool minimal_merge_classes(unsigned k1,
                                    unsigned k2,
                                    std::vector<unsigned> &classpath,
                                    std::vector<unsigned> &cardinalities,
                                    std::vector<unsigned> &queue);

  static void minimal_compress_classpath(std::vector<unsigned> &classpath);

  static std::vector<BlockSystem> non_trivial_transitive(PermGroup const &pg);

  static std::vector<BlockSystem> non_trivial_non_transitive(PermGroup const &pg);

  static std::vector<Block> non_trivial_find_representatives(
    PermSet const &generators,
    std::vector<std::vector<BlockSystem>> const &partial_blocksystems,
    std::vector<unsigned> const &domain_offsets);

  static std::vector<BlockSystem> non_trivial_from_representatives(
    PermSet const &generators,
    std::vector<Block> const &representatives);

  unsigned _degree;
  std::vector<Block> _blocks;
  std::vector<unsigned> _block_indices;
};

std::ostream &operator<<(std::ostream &os, BlockSystem const &bs);

} // namespace mpsym

#endif // _GUARD_BLOCK_SYSTEM_H
