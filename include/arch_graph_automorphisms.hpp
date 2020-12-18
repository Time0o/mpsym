#ifndef GUARD_ARCH_GRAPH_AUTOMORPHISMS_H
#define GUARD_ARCH_GRAPH_AUTOMORPHISMS_H

#include <sstream>
#include <sstream>
#include <string>
#include <vector>

#include "arch_graph_system.hpp"
#include "dump.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "util.hpp"

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

  std::string to_json() const override
  {
    auto bsgs(_automorphisms.bsgs());

    auto sgs(bsgs.strong_generators());
    std::sort(sgs.begin(), sgs.end());

    std::stringstream ss;

    ss << "{\"automorphisms\": ["
       << bsgs.degree() << ","
       << DUMP(bsgs.base()) << ","
       << TRANSFORM_AND_DUMP(std::vector<Perm>(sgs.begin(), sgs.end()),
                             [](Perm const &perm)
                             { return '"' + util::stream(perm) + '"'; })
       << "]}";

    return ss.str();
  }

  unsigned automorphisms_degree() const override
  { return _automorphisms.degree(); }

private:
  PermGroup automorphisms_(AutomorphismOptions const *,
                           internal::timeout::flag) override
  { return _automorphisms; }

  PermGroup _automorphisms;
};

} // namespace internal

} // namespace mpsym

#endif // GUARD_ARCH_GRAPH_AUTOMORPHISMS_H
