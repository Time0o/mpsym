#include <algorithm>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "arch_graph.hpp"
#include "arch_graph_cluster.hpp"
#include "arch_graph_system.hpp"
#include "arch_uniform_super_graph.hpp"
#include "perm.hpp"
#include "perm_group.hpp"
#include "perm_set.hpp"
#include "task_mapping.hpp"
#include "task_orbits.hpp"

namespace mp = boost::multiprecision;
namespace py = pybind11;

using namespace py::literals;

using mpsym::ArchGraph;
using mpsym::ArchGraphCluster;
using mpsym::ArchGraphSystem;
using mpsym::ArchUniformSuperGraph;
using mpsym::TaskMapping;
using mpsym::TaskOrbits;

using mpsym::internal::BSGS;
using mpsym::internal::Perm;
using mpsym::internal::PermGroup;
using mpsym::internal::PermSet;

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

template<typename T = unsigned>
using Set = std::set<T>;

template<typename T>
std::string stream(T const &obj)
{
  std::stringstream ss;
  ss << obj;
  return ss.str();
}

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
    .def_static("from_lua_file", &ArchGraphSystem::from_lua_file,
                "lua_file"_a, "args"_a = std::vector<std::string>())
    .def_static("from_lua", &ArchGraphSystem::from_lua,
                "lua"_a, "args"_a = std::vector<std::string>())
    .def("num_processors", &ArchGraphSystem::num_processors)
    .def("num_channels", &ArchGraphSystem::num_channels)
    .def("automorphisms",
         [](ArchGraphSystem &self)
         { return self.automorphisms(); })
    .def("num_automorphisms",
         [](ArchGraphSystem &self)
         { return self.num_automorphisms(); })
    .def("representative",
         [&](ArchGraphSystem &self, Sequence<> const &mapping)
         {
           return sequence_apply(mapping,
                                 [&](Sequence<> const &mapping)
                                 { return self.repr(mapping); });
         },
         "mapping"_a)
    .def("representative",
         [&](ArchGraphSystem &self,
             Sequence<> const &mapping,
             TaskOrbits *representatives)
         {
           auto num_orbits_old = representatives->num_orbits();

           auto repr(sequence_apply(
             mapping,
             [&](Sequence<> const &mapping)
             { return self.repr(mapping, representatives); }));

           bool repr_is_new = representatives->num_orbits() > num_orbits_old;

           return std::pair<Sequence<>, bool>(repr, repr_is_new);
         },
         "mapping"_a, "representatives"_a)
    .def("orbit",
         [&](ArchGraphSystem &self, Sequence<> const &mapping, bool sorted)
         {
           auto orbit(sequence_multiplex_apply(
             mapping,
             [&](Sequence<> const &mapping)
             { return self.orbit(mapping); }));

           if (sorted)
             std::sort(orbit.begin(), orbit.end());

           return orbit;
         },
         "mapping"_a, "sorted"_a = true);

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
    .def(py::self == py::self)
    .def(py::self != py::self)
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
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def("__len__", &PermGroup::order)
    .def("__iter__",
         [](PermGroup const &self)
         { return py::make_iterator(self.begin(), self.end()); },
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
