#include <algorithm>
#include <iterator>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "arch_graph_automorphisms.h"
#include "arch_graph_cluster.h"
#include "arch_graph_system.h"
#include "arch_uniform_super_graph.h"
#include "dump.h"
#include "perm.h"
#include "perm_set.h"
#include "permlib.h"
#include "profile_parse.h"
#include "profile_util.h"
#include "task_allocation.h"

namespace
{

std::vector<std::string> split_generators(std::string const &gen_str)
{
  auto gen_first = gen_str.find('(');
  auto gen_last = gen_str.rfind(')');

  auto gen_str_trimmed(gen_str.substr(gen_first, gen_last - gen_first + 1u));

  auto res(profile::split(gen_str_trimmed, "),"));

  if (res.size() > 0u) {
    for (auto i = 0u; i < res.size() - 1u; ++i)
      res[i] += ")";
  }

  return res;
}

using gen_type = std::vector<std::vector<std::vector<unsigned>>>;

std::pair<gen_type, unsigned> parse_generators(
  std::vector<std::string> const &gen_strs)
{
  gen_type gens;

  unsigned largest_moved_point = 0u;

  for (auto const &gen_str : gen_strs) {
    gen_type::value_type perm;
    gen_type::value_type::value_type cycle;

    int n_beg = -1;
    for (int i = 0; i < static_cast<int>(gen_str.size()); ++i) {
      char c = gen_str[i];

      switch (c) {
      case '(':
        cycle.clear();
        break;
      case ',':
      case ')':
        {
          int n = profile::stox<int>(gen_str.substr(n_beg, i - n_beg));

          largest_moved_point = std::max(largest_moved_point,
                                         static_cast<unsigned>(n));

          cycle.push_back(n);
          if (c == ')')
            perm.push_back(cycle);

          n_beg = -1;
        }
        break;
      default:
        if (n_beg == -1)
          n_beg = i;
      }
    }

    gens.push_back(perm);
  }

  return {gens, largest_moved_point};
}

cgtl::PermSet convert_generators_mpsym(unsigned degree, gen_type const &gens)
{
  cgtl::PermSet gens_conv;

  for (auto const &gen : gens)
    gens_conv.emplace(degree, gen);

  return gens_conv;
}

permlib::PermSet convert_generators_permlib(unsigned degree,
                                            gen_type const &gens)
{
  std::vector<permlib::Permutation::ptr> gens_conv(gens.size());

  for (auto i = 0u; i < gens.size(); ++i) {
    auto &gen(gens[i]);

    std::stringstream gen_str;
    for (auto j = 0u; j < gen.size(); ++j) {
      auto &cycle(gens[i][j]);
      for (auto k = 0u; k < cycle.size(); ++k) {
        gen_str << cycle[k];

        if (k == cycle.size() - 1) {
          if (j < gen.size() - 1)
            gen_str << ", ";
        } else {
          gen_str << " ";
        }
      }
    }

    gens_conv[i] = permlib::Permutation::ptr(
      new permlib::Permutation(degree, gen_str.str()));
  }

  return {degree, gens_conv};
}

std::tuple<unsigned, unsigned, std::vector<cgtl::TaskAllocation>>
split_task_allocations(std::string const &task_allocations_str,
                       std::string const &regex_str,
                       char delim)
{
  unsigned num_tasks = 0u;

  unsigned min_pe = UINT_MAX;
  unsigned max_pe = 0u;

  std::vector<cgtl::TaskAllocation> task_allocations;

  std::stringstream ss(task_allocations_str);
  std::regex re_task_allocation(regex_str);

  std::string line;
  while (std::getline(ss, line)) {
    std::smatch m;
    if (!std::regex_match(line, m, re_task_allocation))
      throw std::invalid_argument("malformed task allocation expression");

    cgtl::TaskAllocation task_allocation;

    std::string task_allocation_str(m[1]);
    std::size_t pos_begin = 0u;
    std::size_t pos_end;
    unsigned pe;

    for (;;) {
      pos_end = task_allocation_str.find(delim, pos_begin);

      pe = profile::stox<unsigned>(
        pos_end == std::string::npos ?
          task_allocation_str.substr(pos_begin) :
          task_allocation_str.substr(pos_begin, pos_end - pos_begin));

      min_pe = std::min(pe, min_pe);
      max_pe = std::max(pe, max_pe);

      task_allocation.push_back(pe);

      if (pos_end == std::string::npos)
        break;

      pos_begin = pos_end + 1u;
    }

    if (num_tasks == 0u) {
      num_tasks = task_allocation.size();
    } else if (task_allocation.size() != num_tasks) {
      throw std::invalid_argument(
        "currently only equally sized task sets are supported");
    }

    task_allocations.push_back(task_allocation);
  }

  return {min_pe, max_pe, task_allocations};
}

} // anonymous namespace

namespace profile
{

GenericGroup parse_group(std::string const &group_str)
{
  static std::regex re_group("degree:(\\d+),order:(\\d+),gens:(.*)");

  static std::string re_perm = R"((\(\)|(\((\d+,)+\d+\))+))";
  static std::regex re_generators("\\[(" + re_perm + ",)*(" + re_perm + ")?\\]");

  std::smatch m;
  if (!std::regex_match(group_str, m, re_group))
    throw std::invalid_argument("malformed group expression");

  unsigned degree = stox<unsigned>(m[1]);

  unsigned long long order;
  try {
    order = stox<unsigned long long>(m[2]);
  } catch (std::invalid_argument const &) {
    throw std::invalid_argument("group order too large");
  }

  std::string gen_str = m[3];

  if (!std::regex_match(gen_str, re_generators))
    throw std::invalid_argument("malformed generator expression");

  return {degree, order, gen_str};
}

gap::PermSet parse_generators_gap(unsigned degree, std::string const &gen_str)
{ return {degree, gen_str}; }

cgtl::PermSet parse_generators_mpsym(unsigned degree, std::string const &gen_str)
{
   gen_type gen_vect;
   unsigned largest_moved_point;

   std::tie(gen_vect, largest_moved_point) = parse_generators(split_generators(gen_str));

   return convert_generators_mpsym(degree == 0 ? largest_moved_point : degree, gen_vect);
}

permlib::PermSet parse_generators_permlib(unsigned degree, std::string const &gen_str)
{
   gen_type gen_vect;
   unsigned largest_moved_point;

   std::tie(gen_vect, largest_moved_point) = parse_generators(split_generators(gen_str));

   return convert_generators_permlib(degree == 0 ? largest_moved_point : degree, gen_vect);
}

gap::TaskAllocationVector parse_task_allocations_gap(
  std::string const &task_allocations_str)
{
  auto task_allocations(
    split_task_allocations(task_allocations_str, R"((\d+(?: \d+)*))", ' '));

  std::stringstream ss;
  for (auto const &task_allocation : std::get<2>(task_allocations))
    ss << DUMP(task_allocation) << ",\n";

  return {std::get<0>(task_allocations),
          std::get<1>(task_allocations),
          ss.str()};
}

cgtl::TaskAllocationVector parse_task_allocations_mpsym(
  std::string const &task_allocations_str)
{
  auto task_allocations(
    split_task_allocations(task_allocations_str, R"((\d+(?: \d+)*))", ' '));

  return {std::get<0>(task_allocations),
          std::get<1>(task_allocations),
          std::get<2>(task_allocations)};
}

cgtl::TaskAllocationVector parse_task_allocations_gap_to_mpsym(
  std::string const &gap_output_str)
{
  static std::regex re_task_allocations(
    R"(Found \d+ orbit representatives\n((?:.|\n)*))");

  std::smatch m;
  if (!std::regex_search(gap_output_str, m, re_task_allocations))
    throw std::invalid_argument("malformed gap output");

  auto task_allocations(split_task_allocations(
    m[1], R"(.*\[ (\d+(?:, \d+)*) \])", ','));

  return {std::get<0>(task_allocations),
          std::get<1>(task_allocations),
          std::get<2>(task_allocations)};
}

} // namespace profile
