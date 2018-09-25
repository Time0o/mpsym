#ifndef _GUARD_BSGS_H
#define _GUARD_BSGS_H

#include <memory>
#include <ostream>
#include <vector>

#include "perm.h"
#include "schreier_structure.h"

namespace cgtl
{

struct BSGS
{
  std::vector<unsigned> base;
  std::vector<Perm> strong_generators;
  std::vector<std::shared_ptr<SchreierStructure>> schreier_structures;

  bool contains(Perm const &perm) const;
  std::vector<unsigned> orbit(unsigned i) const;
  Perm transversal(unsigned i, unsigned o) const;
  std::vector<Perm> transversals(unsigned i) const;
  std::vector<Perm> stabilizers(unsigned i) const;

  // TODO: add option to keep original generators
  void remove_generators();

  static BSGS solve(std::vector<Perm> const &generators);
};

std::ostream& operator<<(std::ostream& stream, BSGS const &bsgs);

} // namespace cgtl

#endif // _GUARD_BSGS_H
