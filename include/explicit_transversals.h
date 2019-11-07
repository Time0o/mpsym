#ifndef _GUARD_EXPLICIT_TRANSVERSALS_H
#define _GUARD_EXPLICIT_TRANSVERSALS_H

#include <map>
#include <ostream>
#include <vector>

#include "perm.h"
#include "perm_set.h"
#include "schreier_structure.h"

namespace cgtl
{

struct ExplicitTransversals : public SchreierStructure
{
  ExplicitTransversals(unsigned degree)
  : _degree(degree)
  {}

  void create_root(unsigned root) override;
  void create_labels(PermSet const &labels) override;
  void create_edge(
    unsigned origin, unsigned destination, unsigned label) override;

  unsigned root() const override;
  std::vector<unsigned> nodes() const override;
  PermSet labels() const override;

  bool contains(unsigned node) const override;
  bool incoming(unsigned node, Perm const &edge) const override;
  Perm transversal(unsigned origin) const override;

private:
  void dump(std::ostream &os) const override;

  unsigned _degree;
  unsigned _root = 0;
  PermSet _labels;
  std::map<unsigned, Perm> _orbit;
};

} // namespace cgtl

#endif // _GUARD_EXPLICIT_TRANSVERSALS_H
