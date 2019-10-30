#ifndef _GUARD_BSGS_H
#define _GUARD_BSGS_H

#include <cassert>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "perm.h"
#include "schreier_structure.h"

namespace cgtl
{

struct BSGS
{
  BSGS(unsigned degree)
  : degree(degree)
  {}

  unsigned degree;
  std::vector<unsigned> base;
  std::vector<Perm> strong_generators;
  std::vector<std::shared_ptr<SchreierStructure>> schreier_structures;

  std::vector<unsigned> orbit(unsigned i) const;
  Perm transversal(unsigned i, unsigned o) const;
  std::vector<Perm> transversals(unsigned i) const;
  std::vector<Perm> stabilizers(unsigned i) const;

  std::pair<Perm, unsigned> strip(Perm const &perm, unsigned offs = 0) const;
  bool strips_completely(Perm const &perm) const;

  // TODO: allow for other schreier structures
  void extend_base(unsigned bp);
  void update_schreier_structure(unsigned i,
                                 std::vector<Perm> const &strong_generators);
  // TODO: add option to keep original generators
  void remove_generators();
};

std::ostream& operator<<(std::ostream& stream, BSGS const &bsgs);

} // namespace cgtl

#endif // _GUARD_BSGS_H
