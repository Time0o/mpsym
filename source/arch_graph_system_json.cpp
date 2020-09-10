#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "arch_graph.hpp"
#include "arch_graph_automorphisms.hpp"
#include "arch_graph_cluster.hpp"
#include "arch_graph_system.hpp"
#include "arch_uniform_super_graph.hpp"
#include "bsgs.hpp"
#include "dump.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"
#include "util.hpp"

using json = nlohmann::json;

namespace
{

template<typename JSON>
std::shared_ptr<mpsym::ArchGraphSystem>
arch_graph_system_from_json(JSON const &json_)
{
  using mpsym::ArchGraph;
  using mpsym::ArchGraphCluster;
  using mpsym::ArchUniformSuperGraph;

  using mpsym::internal::ArchGraphAutomorphisms;
  using mpsym::internal::BSGS;
  using mpsym::internal::Perm;
  using mpsym::internal::PermGroup;
  using mpsym::internal::PermSet;

  using mpsym::util::parse_perm_set;

  if (!json_.is_object() || json_.size() != 1)
    throw std::logic_error("invalid JSON dictionary");

  auto typestr = json_.begin().key();

  if (typestr == "automorphisms") {
    auto automorphisms(json_["automorphisms"]);

    unsigned degree = automorphisms[0];
    std::vector<unsigned> base = automorphisms[1];
    std::vector<std::string> strong_generators = automorphisms[2];

    PermGroup pg(BSGS(degree, base, parse_perm_set(degree, strong_generators)));

    return std::make_shared<ArchGraphAutomorphisms>(pg);

  } if (typestr == "graph") {
    auto graph(json_["graph"]);

    bool directed = graph[0]["directed"];

    std::vector<std::string> processor_types = graph[1]["processor_types"];
    std::vector<std::string> channel_types = graph[2]["channel_types"];

    using PT = ArchGraph::ProcessorType;
    std::vector<std::pair<PT, std::string>> processors = graph[3]["processors"];

    using CT = std::pair<PT, std::string>;
    std::vector<std::pair<PT, std::vector<CT>>> channels = graph[4]["channels"];

    auto ag(std::make_shared<ArchGraph>(directed));

    for (auto const &pt : processor_types)
      ag->new_processor_type(pt);

    for (auto const &ct : channel_types)
      ag->new_channel_type(ct);

    for (auto const &p : processors)
      ag->add_processor(ag->lookup_processor_type(p.second));

    for (auto const &from : channels) {
      for (auto const &to : from.second)
        ag->add_channel(from.first, to.first, ag->lookup_channel_type(to.second));
    }

    return ag;

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
