#include <algorithm>
#include <chrono>
#include <functional>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <nlohmann/json.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "arch_graph.hpp"
#include "arch_graph_automorphisms.hpp"
#include "arch_graph_cluster.hpp"
#include "arch_graph_system.hpp"
#include "arch_uniform_super_graph.hpp"
#include "nauty_graph.hpp"
#include "parse.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"
#include "task_mapping.hpp"
#include "task_orbits.hpp"
#include "timeout.hpp"
#include "util.hpp"

namespace mp = boost::multiprecision;

namespace py = pybind11;

using namespace py::literals;

using json = nlohmann::json;

using mpsym::ArchGraph;
using mpsym::ArchGraphCluster;
using mpsym::ArchGraphSystem;
using mpsym::ArchUniformSuperGraph;
using mpsym::TaskMapping;
using mpsym::TaskOrbits;

using mpsym::internal::ArchGraphAutomorphisms;
using mpsym::internal::BSGS;
using mpsym::internal::NautyGraph;
using mpsym::internal::Perm;
using mpsym::internal::PermGroup;
using mpsym::internal::PermSet;

using mpsym::util::parse_perm;
using mpsym::util::stream;

using mpsym::internal::timeout::flag;
using mpsym::internal::timeout::run_abortable_with_timeout;
using mpsym::internal::timeout::TimeoutError;

namespace
{

template<typename T = unsigned>
using Sequence = std::vector<T>;

template<typename T>
using contained_type =
  typename std::remove_reference<decltype(*std::declval<T>().begin())>::type;

template<typename T>
Sequence<contained_type<T>> to_sequence(T const &obj)
{ return Sequence<contained_type<T>>(obj.begin(), obj.end()); }

template<typename T = unsigned>
Sequence<T> to_sequence(py::iterator it)
{
  Sequence<T> seq;

  while (it != py::iterator::sentinel()) {
    try {
      seq.push_back(it->cast<T>());
    } catch (py::cast_error const &) {
      throw std::runtime_error("failed to convert iterable to matching sequence");
    }

    ++it;
  }

  return seq;
}

template<typename T = unsigned>
Sequence<T> inc_sequence(Sequence<T> seq)
{
  std::transform(seq.begin(), seq.end(), seq.begin(),
                 [](T x){ return x + 1; });

  return seq;
}

template<typename T = unsigned>
Sequence<T> dec_sequence(Sequence<T> seq)
{
  std::transform(seq.begin(), seq.end(), seq.begin(),
                 [](T x){ return x - 1; });

  return seq;
}

template<typename T, typename FUNC>
Sequence<T> sequence_apply(Sequence<T> const &seq, FUNC &&f)
{ return dec_sequence(to_sequence(f(inc_sequence(seq)))); }

template<typename T, typename FUNC>
Sequence<Sequence<T>> sequence_multiplex_apply(Sequence<T> const &seq, FUNC &&f)
{
  Sequence<Sequence<T>> seqs;
  for (auto const &seq_ : f(inc_sequence(seq)))
    seqs.push_back(dec_sequence(to_sequence(seq_)));

  return seqs;
}

template<typename U, typename T, typename FUNC>
Sequence<U> transform_sequence(Sequence<T> const &seq, FUNC &&f)
{
  Sequence<U> seq_transformed(seq.size());
  std::transform(seq.begin(), seq.end(), seq_transformed.begin(), f);
  return seq_transformed;
}

template<typename T>
py::tuple sequence_to_tuple(Sequence<T> const &seq)
{
  py::tuple t = py::cast(seq);
  return t;
}

template<typename T = unsigned>
using Set = std::set<T>;

template<typename T>
T arch_graph_json_cast(ArchGraphSystem const &self, std::string const &key)
{
  auto j(json::parse(self.to_json()));
  T ret = j["graph"][key];
  return ret;
}

#define AGS_TIMEOUT(what, method) \
  run_abortable_with_timeout( \
    what, \
    std::chrono::duration<double>(timeout), \
    [&](flag aborted) \
    { return self.method(nullptr, aborted); })

#define AGS_TIMEOUT_WITH_ARGS(what, method, ...) \
  run_abortable_with_timeout( \
    what, \
    std::chrono::duration<double>(timeout), \
    [&](flag aborted) \
    { return self.method(__VA_ARGS__, nullptr, aborted); })

} // anonymous namespace

namespace pybind11
{

namespace detail
{

template <>
struct type_caster<mp::cpp_int> {
  PYBIND11_TYPE_CASTER(mp::cpp_int, _("cpp_int"));

  static handle cast(mp::cpp_int const &src, return_value_policy, handle)
  { return PyLong_FromString(stream(src).c_str(), nullptr, 10); }
};

} // namespace detail

} // namespace pybind11

#define PYBIND11_MODULE_(name, m) PYBIND11_MODULE(name, m)

PYBIND11_MODULE_(PYTHON_MODULE, m)
{
  m.doc() = DESCRIPTION;
  m.attr("__version__") = VERSION;

  // ArchGraphSystem
  py::class_<ArchGraphSystem,
             std::shared_ptr<ArchGraphSystem>>(m, "ArchGraphSystem")
    .def("__repr__", &ArchGraphSystem::to_json)
    .def_static("from_lua", &ArchGraphSystem::from_lua,
                "lua"_a, "args"_a = std::vector<std::string>())
    .def_static("from_lua_file", &ArchGraphSystem::from_lua_file,
                "lua_file"_a, "args"_a = std::vector<std::string>())
    .def_static("from_nauty",
                [](int vertices,
                   std::map<int, std::vector<int>> const &adjacencies,
                   int vertices_reduced,
                   bool directed,
                   std::vector<std::set<int>> const &coloring)
                {
                  // validate number of vertices
                  if (vertices <= 0)
                    throw std::logic_error("number of vertices must be non-negative");

                  // validate adjacencies
                  std::set<std::pair<int, int>> edges;
                  std::set<std::pair<int, int>> edge_parities;

                  for (auto const &p : adjacencies) {
                    int from = p.first;

                    if (from >= vertices)
                      throw std::logic_error("vertex index out of range");

                    for (int to : p.second) {
                      if (to >= vertices)
                        throw std::logic_error("vertex index out of range");

                      auto edge(std::make_pair(from, to));
                      auto reverse_edge(std::make_pair(to, from));

                      if (!edges.insert(edge).second)
                        throw std::logic_error("duplicate edges");
                    }
                  }

                  // validate coloring
                  if (!coloring.empty()) {
                    std::set<int> tmp;
                    for (auto const &p : coloring)
                      tmp.insert(p.begin(), p.end());

                    if (static_cast<int>(tmp.size()) != vertices
                        || *tmp.begin() != 0
                        || *tmp.rbegin() != vertices - 1) {

                      throw std::logic_error("invalid coloring");
                    }
                  }

                  // construct graph
                  if (vertices_reduced == 0)
                    vertices_reduced = vertices;

                  NautyGraph g{vertices, directed};

                  g.add_edges(adjacencies);

                  if (!coloring.empty()) {
                    g.set_partition(
                      transform_sequence<std::vector<int>>(
                        coloring,
                        [](std::set<int> const &p)
                        { return std::vector<int>(p.begin(), p.end()); }));
                  }

                  // extract automorphisms
                  return std::make_shared<ArchGraphAutomorphisms>(
                    PermGroup(vertices_reduced, g.automorphism_generators()));
                },
                "vertices"_a,
                "adjacencies"_a,
                "vertices_reduced"_a = 0,
                "directed"_a = true,
                "coloring"_a = std::vector<int>())
    .def_static("from_json", &ArchGraphSystem::from_json,
                "json"_a)
    .def_static("from_json_file", &ArchGraphSystem::from_json_file,
                "json_file"_a)
    .def("to_json", &ArchGraphSystem::to_json)
    .def("processor_types",
         [](ArchGraphSystem const &self)
         {
           using T = std::vector<std::string>;
           return arch_graph_json_cast<T>(self, "processor_types");
         })
    .def("channel_types",
         [](ArchGraphSystem const &self)
         {
           using T = std::vector<std::string>;
           return arch_graph_json_cast<T>(self, "channel_types");
         })
    .def("processors",
         [](ArchGraphSystem const &self)
         {
           using T = std::pair<unsigned, std::string>;
           using U = std::vector<T>;
           return arch_graph_json_cast<U>(self, "processors");
         })
    .def("channels",
         [](ArchGraphSystem const &self)
         {
           using T = std::pair<unsigned, std::string>;
           using U = std::map<unsigned, std::vector<T>>;
           return arch_graph_json_cast<U>(self, "channels");
         })
    .def("num_processors", &ArchGraphSystem::num_processors)
    .def("num_channels", &ArchGraphSystem::num_channels)
    .def("initialize",
         [&](ArchGraphSystem &self, double timeout)
         { return AGS_TIMEOUT("initialize", init_repr); },
         "timeout"_a = 0.0)
    .def("num_automorphisms",
         [](ArchGraphSystem &self, double timeout)
         { return AGS_TIMEOUT("num_automorphisms", num_automorphisms); },
         "timeout"_a = 0.0)
    .def("automorphisms",
         [](ArchGraphSystem &self, double timeout)
         { return AGS_TIMEOUT("automorphisms", automorphisms); },
         "timeout"_a = 0.0)
    .def("expand_automorphisms", &ArchGraphSystem::expand_automorphisms)
    .def("representative",
         [&](ArchGraphSystem &self, Sequence<> const &mapping, double timeout)
         {
           auto repr(AGS_TIMEOUT_WITH_ARGS("representative",
                                           repr,
                                           inc_sequence(mapping)));

           return sequence_to_tuple(dec_sequence(to_sequence(repr)));
         },
         "mapping"_a, "timeout"_a = 0.0)
    .def("representative",
         [&](ArchGraphSystem &self,
             Sequence<> const &mapping,
             TaskOrbits &representatives,
             double timeout)
         {
           TaskMapping repr;
           bool orbit_new;
           unsigned orbit_index;

           std::tie(repr, orbit_new, orbit_index) =
             AGS_TIMEOUT_WITH_ARGS("representative",
                                   repr,
                                   inc_sequence(mapping),
                                   representatives);

           return std::make_tuple(
             sequence_to_tuple(dec_sequence(to_sequence(repr))),
             orbit_new,
             orbit_index);
         },
         "mapping"_a, "representatives"_a, "timeout"_a = 0.0)
    .def("orbit",
         [&](ArchGraphSystem &self,
             Sequence<> const &mapping,
             bool sorted,
             double timeout)
         {
           auto orb(sequence_multiplex_apply(
             mapping,
             [&](Sequence<> const &mapping)
             {
               return run_abortable_with_timeout(
                 "orbit",
                 std::chrono::duration<double>(timeout),
                 [&](flag aborted)
                 { return self.orbit(mapping, aborted); });
             }));

           if (sorted)
             std::sort(orb.begin(), orb.end());

           return transform_sequence<py::tuple>(
             orb,
             [](Sequence<> const &mapping)
             { return sequence_to_tuple(mapping); });
         },
         "mapping"_a, "sorted"_a = true, "timeout"_a = 0.0);

  // ArchGraphAutomorphisms
  py::class_<ArchGraphAutomorphisms,
             ArchGraphSystem,
             std::shared_ptr<ArchGraphAutomorphisms>>(m, "ArchGraphAutomorphisms")
    .def(py::init<PermGroup>(), "automorphisms"_a)
    .def(py::pickle(
        [](ArchGraphAutomorphisms &self)
        { return self.to_json(); },
        [](std::string const &json)
        {
          return std::dynamic_pointer_cast<ArchGraphAutomorphisms>(
            ArchGraphSystem::from_json(json));
        }));

  // ArchGraph
  py::class_<ArchGraph,
             ArchGraphSystem,
             std::shared_ptr<ArchGraph>>(m, "ArchGraph")
    .def(py::init<bool>(), "directed"_a = true)
    .def(py::pickle(
        [](ArchGraph &self)
        { return self.to_json(); },
        [](std::string const &json)
        {
          return std::dynamic_pointer_cast<ArchGraph>(
            ArchGraphSystem::from_json(json));
        }))
    .def("directed", &ArchGraph::directed)
    .def("add_processor",
         (unsigned(ArchGraph::*)(std::string const &))
         &ArchGraph::add_processor,
         "pl"_a)
    .def("add_processors",
         (unsigned(ArchGraph::*)(unsigned, std::string const &))
         &ArchGraph::add_processors,
         "n"_a, "pl"_a)
    .def("add_channel",
         (void(ArchGraph::*)(unsigned, unsigned, std::string const &))
         &ArchGraph::add_channel,
         "pe1"_a, "pe2"_a, "cl"_a)
    .def("add_channels",
         (void(ArchGraph::*)(ArchGraph::ChannelDict<std::string> const &))
         &ArchGraph::add_channels<std::string>,
         "channels"_a)
    .def("add_channels",
         (void(ArchGraph::*)(ArchGraph::ChannelDict<> const &, std::string const &))
         &ArchGraph::add_channels<std::string const &>,
         "channels"_a, "ct"_a)
    .def("fully_connect",
         (void(ArchGraph::*)(std::string const &))
         &ArchGraph::fully_connect,
         "cl"_a)
    .def("fully_connect",
         (void(ArchGraph::*)(std::string const &, std::string const &))
         &ArchGraph::fully_connect,
         "pl"_a, "cl"_a)
    .def("fully_connect",
         [](ArchGraph &self, py::iterable it, std::string const &cl)
         {
           auto processors(to_sequence<>(py::make_iterator(it)));

           for (unsigned pe1 : processors) {
             for (unsigned pe2 : processors) {
               if (pe1 != pe2)
                 self.add_channel(pe1, pe2, cl);
             }
           }
         },
         "processors"_a, "cl"_a)
    .def("self_connect",
         (void(ArchGraph::*)(std::string const &))
         &ArchGraph::self_connect,
         "cl"_a)
    .def("self_connect",
         (void(ArchGraph::*)(std::string const &, std::string const &))
         &ArchGraph::self_connect,
         "pl"_a, "cl"_a)
    .def("self_connect",
         [](ArchGraph &self, py::iterator it, std::string const &cl)
         {
           auto processors(to_sequence<>(it));

           for (unsigned pe : processors)
             self.add_channel(pe, pe, cl);
         },
         "processors"_a, "cl"_a);

  // ArchGraphCluster
  py::class_<ArchGraphCluster,
             ArchGraphSystem,
             std::shared_ptr<ArchGraphCluster>>(m, "ArchGraphCluster")
    .def(py::init<>())
    .def(py::pickle(
        [](ArchGraphCluster &self)
        { return self.to_json(); },
        [](std::string const &json)
        {
          return std::dynamic_pointer_cast<ArchGraphCluster>(
            ArchGraphSystem::from_json(json));
        }))
    .def("add_subsystem", &ArchGraphCluster::add_subsystem, "subsystem"_a)
    .def("num_subsystems", &ArchGraphCluster::num_subsystems);

  // ArchUniformSuperGraph
  py::class_<ArchUniformSuperGraph,
             ArchGraphSystem,
             std::shared_ptr<ArchUniformSuperGraph>>(m, "ArchUniformSuperGraph")
    .def(py::init<std::shared_ptr<ArchGraphSystem>,
                  std::shared_ptr<ArchGraphSystem>>(),
         "super_graph"_a, "proto"_a)
    .def(py::pickle(
        [](ArchUniformSuperGraph &self)
        { return self.to_json(); },
        [](std::string const &json)
        {
          return std::dynamic_pointer_cast<ArchUniformSuperGraph>(
            ArchGraphSystem::from_json(json));
        }));

  // TaskOrbits
  py::class_<TaskOrbits>(m, "Representatives")
    .def(py::init<>())
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def("__len__", &TaskOrbits::num_orbits)
    .def("__iter__",
         [](TaskOrbits const &orbits)
         { return py::make_iterator(orbits.begin(), orbits.end()); })
    .def("__contains__",
         [](TaskOrbits const &orbits, Sequence<> const &mapping)
         { return orbits.is_repr(inc_sequence(mapping)); },
         "mapping"_a);

  // Perm
  py::class_<Perm>(m, "Perm")
    .def(py::init<unsigned>(), "degree"_a = 1)
    .def(py::init(
           [](Sequence<> const &v)
           {
             if (v.empty())
                throw std::logic_error("invalid permutation");

             auto max = *std::max_element(v.begin(), v.end());
             if (v.size() != max + 1)
               throw std::logic_error("invalid permutation");

             Set<> s(v.begin(), v.end());
             if (s.size() != max + 1 || *s.begin() != 0)
               throw std::logic_error("invalid permutation");

             return Perm(inc_sequence(v));
           }),
         "perm"_a)
    .def(py::init(
           [](unsigned degree, Sequence<Sequence<>> const &cycles)
           {
             std::unordered_set<unsigned> cycles_flattened;

             for (auto const &cycle : cycles) {
                for (unsigned x : cycle) {
                  if (x > degree || !cycles_flattened.insert(x).second)
                    throw std::invalid_argument("invalid permutation");
                }
             }

             return Perm(degree, cycles);
           }),
           "degree"_a, "cycles"_a)
    .def(py::init(
           [](unsigned degree, std::string cycles)
           {
             static std::regex re_perm(R"((\(\)|(\(( *\d+,)+ *\d+ *\))+))");

             cycles.erase(std::remove(cycles.begin(), cycles.end(), ' ' ), cycles.end());

             if (!std::regex_match(cycles, re_perm))
               throw std::invalid_argument("invalid permutation string");

             return parse_perm(degree, cycles);
           }),
           "degree"_a, "cycles"_a)
    .def("__eq__",
         [](Perm const &self, Perm const &other)
         {
           if (self.degree() != other.degree())
             throw std::invalid_argument("can only compare permutations of equal degree");

           return self == other;
         })
    .def("__ne__",
         [](Perm const &self, Perm const &other)
         {
           if (self.degree() != other.degree())
             throw std::invalid_argument("can only compare permutations of equal degree");

           return self != other;
         })
    .def(hash(py::self))
    .def("__getitem__",
         [](Perm const &self, unsigned x)
         {
           if (x > self.degree() - 1)
             throw std::out_of_range("not in domain");

           return self[x + 1] - 1;
         },
         "x"_a)
    .def("__invert__", &Perm::operator~)
    .def("__mul__",
         [](Perm const &self, Perm const &other)
         {
           if (self.degree() != other.degree())
             throw std::logic_error("permutation degrees do not match");

           return self * other;
         },
         "other"_a)
    .def("__rmul__",
         [](Perm const &self, Perm const &other)
         {
           if (self.degree() != other.degree())
             throw std::logic_error("permutation degrees do not match");

           return other * self;
         },
         "other"_a)
    .def("__bool__",
         [](Perm const &self)
         { return !self.id(); })
    .def("__repr__",
         [](Perm const &self)
         { return stream(self); })
    .def("degree", &Perm::degree);

  py::implicitly_convertible<Sequence<>, Perm>();

  // PermGroup
  py::class_<PermGroup>(m, "PermGroup")
    .def(py::init<unsigned>(), "degree"_a = 1)
    .def(py::init(
           [](Sequence<Perm> const &generators_)
           {
             if (generators_.empty())
               throw std::logic_error("generating set must not be empty");

             unsigned degree = generators_[0].degree();
             for (std::size_t i = 1; i < generators_.size(); ++i) {
                if (generators_[i].degree() != degree)
                  throw std::logic_error("mismatched generator degrees");
             }

             PermSet generators(generators_.begin(), generators_.end());

             return PermGroup(degree, generators);
           }),
         "generators"_a)
    .def_static("symmetric", PermGroup::symmetric, "degree"_a)
    .def_static("cyclic", PermGroup::cyclic, "degree"_a)
    .def_static("dihedral", PermGroup::dihedral, "degree"_a)
    .def_static("direct_product",
                [](py::iterable it)
                {
                  auto groups(to_sequence<PermGroup>(py::make_iterator(it)));

                  return PermGroup::direct_product(groups.begin(), groups.end());
                },
                "groups"_a)
    .def_static("wreath_product",
                [](PermGroup const &lhs, PermGroup const &rhs)
                { return PermGroup::wreath_product(lhs, rhs); },
                "lhs"_a, "rhs"_a)
    .def("__eq__",
         [](PermGroup const &self, PermGroup const &other)
         {
           if (self.degree() != other.degree())
             throw std::invalid_argument("can only compare permutation groups of equal degree");

           return self == other;
         })
    .def("__ne__",
         [](PermGroup const &self, PermGroup const &other)
         {
           if (self.degree() != other.degree())
             throw std::invalid_argument("can only compare permutation groups of equal degree");

           return self != other;
         })
    .def("__len__", &PermGroup::order)
    .def("__iter__",
         [](PermGroup const &self)
         { return py::make_iterator<py::return_value_policy::copy>(self.begin(), self.end()); },
         py::keep_alive<0, 1>())
    .def("__contains__",
         [](PermGroup const &self, Perm const &p)
         {
           if (p.degree() != self.degree())
             throw std::logic_error("mismatched degrees");

           return self.contains_element(p);
         },
         "perm"_a)
    .def("__bool__",
         [](PermGroup const &self)
         { return !self.is_trivial(); })
    .def("__repr__",
         [](PermGroup const &self)
         { return stream(self.generators()); })
    .def("degree", &PermGroup::degree)
    .def("bsgs",
         [&](PermGroup const &self)
         {
           auto bsgs(self.bsgs());

           return std::make_pair<BSGS::Base, Sequence<Perm>>(
             bsgs.base(),
             to_sequence(bsgs.strong_generators()));
         })
    .def("generators",
         [&](PermGroup const &self, bool sorted)
         {
           auto generators(self.generators());

           if (sorted)
             std::sort(generators.begin(), generators.end());

           return to_sequence(generators);
         },
         "sorted"_a = true)
    .def("is_symmetric", &PermGroup::is_symmetric)
    .def("is_alternating", &PermGroup::is_alternating)
    .def("is_transitive", &PermGroup::is_transitive);

  py::implicitly_convertible<Sequence<Perm>, PermGroup>();
}
