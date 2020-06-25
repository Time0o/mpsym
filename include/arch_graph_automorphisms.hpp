#ifndef GUARD_ARCH_GRAPH_AUTOMORPHISMS_H
#define GUARD_ARCH_GRAPH_AUTOMORPHISMS_H

#include <sstream>
#include <string>

#include "arch_graph_system.hpp"
#include "bsgs.hpp"
#include "perm_group.hpp"

namespace mpsym
{

namespace internal
{

class ArchGraphAutomorphisms : public ArchGraphSystem
{
public:
  ArchGraphAutomorphisms(PermGroup const &automorphisms)
  : _automorphisms(automorphisms)
  {}

  virtual ~ArchGraphAutomorphisms() = default;

  std::string to_gap() const override
  {
    std::stringstream ss;
    ss << _automorphisms.generators();

    auto generators_str(ss.str());
    generators_str.front() = '[';
    generators_str.back() = ']';

    return "Group(" + generators_str + ")";
  }

private:
  PermGroup automorphisms_(AutomorphismOptions const *) override
  { return _automorphisms; }

  PermGroup _automorphisms;
};

} // namespace internal

} // namespace mpsym

#endif // GUARD_ARCH_GRAPH_AUTOMORPHISMS_H
