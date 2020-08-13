#include <memory>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "arch_graph.hpp"
#include "arch_graph_cluster.hpp"
#include "arch_graph_system.hpp"
#include "arch_uniform_super_graph.hpp"
#include "task_mapping.hpp"

#define PYBIND11_MODULE_(name, m) PYBIND11_MODULE(name, m)

PYBIND11_MODULE_(PYMPSYM, m) {
  namespace py = pybind11;

  using namespace py::literals;
  using namespace std::placeholders;

  using mpsym::ArchGraph;
  using mpsym::ArchGraphCluster;
  using mpsym::ArchGraphSystem;
  using mpsym::ArchUniformSuperGraph;
  using mpsym::TaskMapping;

  m.doc() = DESCRIPTION;
  m.attr("__version__") = VERSION;

  // ArchGraphSystem
  py::class_<ArchGraphSystem,
             std::shared_ptr<ArchGraphSystem>>(m, "ArchGraphSystem")
    .def_static("from_lua_file", &ArchGraphSystem::from_lua_file,
                "lua_file"_a, "args"_a = std::vector<std::string>())
    .def_static("from_lua", &ArchGraphSystem::from_lua,
                "lua"_a, "args"_a = std::vector<std::string>())
    .def("num_processors", &ArchGraphSystem::num_processors)
    .def("num_channels", &ArchGraphSystem::num_channels)
    .def("representative",
         [](ArchGraphSystem &self, std::vector<unsigned> const &mapping)
         {
            auto repr(self.repr(mapping));
            return std::vector<unsigned>(repr.begin(), repr.end());
         },
         "mapping"_a);

  // ArchGraph
  py::class_<ArchGraph,
             ArchGraphSystem,
             std::shared_ptr<ArchGraph>>(m, "ArchGraph")
    .def(py::init<>())
    .def("new_processor_type", &ArchGraph::new_processor_type, "pl"_a = "")
    .def("new_channel_type", &ArchGraph::new_channel_type, "cl"_a = "")
    .def("add_processor", &ArchGraph::add_processor, "pe"_a)
    .def("add_channel", &ArchGraph::add_channel, "pe1"_a, "pe2"_a, "ch"_a);

  // ArchGraphCluster
  py::class_<ArchGraphCluster,
             ArchGraphSystem,
             std::shared_ptr<ArchGraphCluster>>(m, "ArchGraphCluster")
    .def(py::init<>())
    .def("add_subsystem", &ArchGraphCluster::add_subsystem, "subsystem"_a)
    .def("num_subsystems", &ArchGraphCluster::num_subsystems);

  // ArchUniformSuperGraph
  py::class_<ArchUniformSuperGraph,
             ArchGraphSystem,
             std::shared_ptr<ArchUniformSuperGraph>>(m, "ArchUniformSuperGraph")
    .def(py::init<std::shared_ptr<ArchGraphSystem>,
                  std::shared_ptr<ArchGraphSystem>>(),
         "super_graph"_a, "proto"_a);
}
