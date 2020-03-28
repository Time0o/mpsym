#ifndef _GUARD_SCHREIER_TREE_H
#define _GUARD_SCHREIER_TREE_H

#include <map>
#include <ostream>
#include <vector>

#include "perm.h"
#include "perm_set.h"
#include "schreier_structure.h"

namespace mpsym
{

struct SchreierTree : public SchreierStructure
{
  SchreierTree(unsigned degree, unsigned root, PermSet const &labels)
  : _degree(degree),
    _root(root),
    _labels(labels)
  {}

  void add_label(Perm const &label)
  { _labels.insert(label); }

  void create_edge(unsigned origin,
                   unsigned destination,
                   unsigned label) override;

  unsigned root() const override;
  std::vector<unsigned> nodes() const override;
  PermSet labels() const override;

  bool contains(unsigned node) const override;
  bool incoming(unsigned node, Perm const &edge) const override;
  Perm transversal(unsigned origin) const override;

private:
  void dump(std::ostream &os) const override;

  unsigned _degree;
  unsigned _root;
  std::map<unsigned, unsigned> _edges;
  PermSet _labels;
  std::map<unsigned, unsigned> _edge_labels;
};

} // namespace mpsym

#endif // _GUARD_SCHREIER_TREE_H
