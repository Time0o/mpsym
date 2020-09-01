#include <sstream>
#include <string>

#include "arch_graph_system.hpp"
#include "dump.hpp"
#include "perm.hpp"

namespace mpsym
{

std::string
ArchGraphSystem::to_json()
{
  using mpsym::internal::Perm;

  auto bsgs(automorphisms().bsgs());

  std::stringstream ss;

  ss << "{\"automorphisms\": ["
     << bsgs.degree() << ","
     << DUMP(bsgs.base()) << ","
     << TRANSFORM_AND_DUMP(bsgs.strong_generators(),
                           [](Perm const &perm){ return perm.vect(); }) << "]}";

  return ss.str();
}

} // namespace mpsym
