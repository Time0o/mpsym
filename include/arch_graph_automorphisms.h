#ifndef _GUARD_ARCH_GRAPH_AUTOMORPHISMS_H
#define _GUARD_ARCH_GRAPH_AUTOMORPHISMS_H

#include "arch_graph_system.h"
#include "perm_group.h"

namespace cgtl
{

class ArchGraphAutomorphisms : public ArchGraphSystem
{
public:
  ArchGraphAutomorphisms(PermGroup const &automorphisms)
  : _automorphisms(automorphisms)
  {}

private:
  PermGroup update_automorphisms() override
  { return _automorphisms; }

  PermGroup _automorphisms;
};

} // namespace cgtl

#endif // _GUARD_ARCH_GRAPH_AUTOMORPHISMS_H
