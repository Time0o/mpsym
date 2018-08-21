#ifndef _GUARD_BLOCK_SYSTEM_H
#define _GUARD_BLOCK_SYSTEM_H

#include <ostream>
#include <vector>

#include "perm.h"

namespace cgtl
{

class BlockSystem
{
friend std::ostream& operator<<(std::ostream& stream, BlockSystem const &bs);

public:
  typedef std::vector<std::vector<unsigned>>::const_iterator const_iterator;

  BlockSystem(std::vector<unsigned> classes);

  unsigned degree() const { return _n; }
  unsigned size() const { return static_cast<unsigned>(_blocks.size()); }

  std::vector<unsigned> const& operator[](unsigned const i) const;
  const_iterator begin() const;
  const_iterator end() const;

  static BlockSystem minimal(std::vector<Perm> const &generators,
                             std::vector<unsigned> const &initial_class);

private:
  unsigned _n;
  std::vector<std::vector<unsigned>> _blocks;
};

std::ostream& operator<<(std::ostream& stream, BlockSystem const &bs);

} // namespace cgtl

#endif // _GUARD_BLOCK_SYSTEM_H
