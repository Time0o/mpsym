#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "arch_graph_automorphisms.hpp"
#include "arch_graph_cluster.hpp"
#include "arch_graph_system.hpp"
#include "arch_uniform_super_graph.hpp"
#include "bsgs.hpp"
#include "dump.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"

using json = nlohmann::json;

namespace
{

template<typename JSON>
std::shared_ptr<mpsym::ArchGraphSystem>
arch_graph_system_from_json(JSON const &json_)
{
  using mpsym::ArchGraphCluster;
  using mpsym::ArchUniformSuperGraph;

  using mpsym::internal::ArchGraphAutomorphisms;
  using mpsym::internal::BSGS;
  using mpsym::internal::Perm;
  using mpsym::internal::PermGroup;
  using mpsym::internal::PermSet;

  if (!json_.is_object() || json_.size() != 1)
    throw std::logic_error("invalid JSON dictionary");

  auto typestr = json_.begin().key();

  if (typestr == "automorphisms") {
    auto automorphisms(json_["automorphisms"]);

    unsigned degree = automorphisms[0];
    std::vector<unsigned> base = automorphisms[1];
    std::vector<std::vector<unsigned>> strong_generators_ = automorphisms[2];

    PermSet strong_generators;
    for (auto const &vect : strong_generators_)
      strong_generators.emplace(vect);

    PermGroup pg(BSGS(degree, base, strong_generators));

    return std::make_shared<ArchGraphAutomorphisms>(pg);

  } else if (typestr == "cluster") {
    auto cluster(json_["cluster"]);

    auto agc(std::make_shared<ArchGraphCluster>());
    for (auto const &ags : cluster)
      agc->add_subsystem(arch_graph_system_from_json(ags));

    return agc;

  } else if (typestr == "super_graph") {
    auto super_graph_(json_["super_graph"]);

    auto proto(arch_graph_system_from_json(super_graph_[0]));
    auto super_graph(arch_graph_system_from_json(super_graph_[1]));

    return std::make_shared<ArchUniformSuperGraph>(super_graph, proto);

  } else {
    throw std::logic_error("invalid JSON dictionary");
  }
}

} // anonymous namespace

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

std::shared_ptr<ArchGraphSystem>
ArchGraphSystem::from_json(std::string const &json_)
{
  try {
    return arch_graph_system_from_json(json::parse(json_));
  } catch (std::runtime_error const &e) {
    throw std::runtime_error("failed to parse JSON: " + std::string(e.what()));
  }
}

} // namespace mpsym
