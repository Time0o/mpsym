#include <algorithm>
#include <regex>
#include <string>
#include <vector>
#include <sstream>
#include <tuple>

#include "dump.h"
#include "perm.h"
#include "perm_set.h"
#include "permlib.h"
#include "profile_parse.h"
#include "profile_utility.h"
#include "task_mapping.h"

namespace
{

using gen_type = std::vector<std::vector<std::vector<unsigned>>>;

std::vector<std::string> split_generators(std::string const &gen_str)
{
  std::vector<std::string> gen_strs;

  std::size_t pos = 1, nextpos;
  for (;;) {
    nextpos = gen_str.find("),", pos);

    if (nextpos == std::string::npos) {
      gen_strs.push_back(gen_str.substr(pos, gen_str.size() - pos - 1));
      break;
    } else {
      gen_strs.push_back(gen_str.substr(pos, nextpos - pos + 1));
    }

    pos = nextpos + 2;
  }

  return gen_strs;
}

std::pair<unsigned, gen_type> parse_generators(
  std::vector<std::string> const &gen_strs)
{
  unsigned degree = 0;

  gen_type gens;

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
          int n = stox<int>(gen_str.substr(n_beg, i - n_beg));

          degree = std::max(degree, static_cast<unsigned>(n));

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

  return {degree, gens};
}

cgtl::PermSet convert_generators_mpsym(unsigned degree, gen_type const &gens)
{
  cgtl::PermSet gens_conv;

  for (auto const &gen : gens)
    gens_conv.emplace(degree, gen);

  return gens_conv;
}

std::vector<permlib::Permutation::ptr> convert_generators_permlib(
  unsigned degree, gen_type const &gens)
{
  (void)degree;

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

  return gens_conv;
}

std::vector<std::vector<unsigned>> split_task_allocations(
  std::string const &task_allocations_str)
{
  static std::regex re_task_allocation(R"(\d+( \d+)*)");

  unsigned num_tasks = 0u;
  std::vector<std::vector<unsigned>> task_allocations;

  std::stringstream ss(task_allocations_str);

  std::string line;
  while (std::getline(ss, line)) {
    if (!std::regex_match(line, re_task_allocation))
      throw std::invalid_argument("malformed task allocation expression");

    std::vector<unsigned> task_allocation;

    std::size_t pos_begin = 0u;
    std::size_t pos_end;
    unsigned pe;

    while ((pos_end = line.find(' ', pos_begin)) != std::string::npos) {
      pe = stox<unsigned>(line.substr(pos_begin, pos_end - pos_begin));
      task_allocation.push_back(pe);

      pos_begin = pos_end + 1u;
    }

    pe = stox<unsigned>(line.substr(pos_begin));
    task_allocation.push_back(pe);

    if (num_tasks == 0u) {
      num_tasks = task_allocation.size();
    } else if (task_allocation.size() != num_tasks) {
      throw std::invalid_argument(
        "currently only equally sized task sets are supported");
    }

    task_allocations.push_back(task_allocation);
  }

  return task_allocations;
}

} // namespace

std::tuple<unsigned, unsigned, std::string> parse_group(
  std::string const &group_str)
{
  static std::regex re_group("degree:(\\d+),order:(\\d+),gens:(.*)");

  static std::string re_perm = R"((\(\)|(\((\d+,)+\d+\))+))";
  static std::regex re_generators("\\[(" + re_perm + ",)*(" + re_perm + ")?\\]");

  std::smatch m;
  if (!std::regex_match(group_str, m, re_group))
    throw std::invalid_argument("malformed group expression");

  unsigned degree = stox<unsigned>(m[1]);
  unsigned order = stox<unsigned>(m[2]);
  std::string gen_str = m[3];

  if (!std::regex_match(gen_str, re_generators))
    throw std::invalid_argument("malformed generator expression");

  return std::make_tuple(degree, order, gen_str);
}

std::string parse_generators_gap(std::string const &gen_str)
{ return gen_str; }

cgtl::PermSet parse_generators_mpsym(std::string const &gen_str)
{
   unsigned degree;
   gen_type gen_vect;
   std::tie(degree, gen_vect) = parse_generators(split_generators(gen_str));

   return convert_generators_mpsym(degree, gen_vect);
}

std::vector<permlib::Permutation::ptr> parse_generators_permlib(
  std::string const &gen_str)
{
   unsigned degree;
   gen_type gen_vect;
   std::tie(degree, gen_vect) = parse_generators(split_generators(gen_str));

   return convert_generators_permlib(degree, gen_vect);
}

std::string parse_task_allocations_gap(std::string const &task_allocations_str)
{
  auto task_allocations(split_task_allocations(task_allocations_str));

  std::stringstream ss;
  for (auto const &task_allocation : task_allocations)
    ss << dump::dump(task_allocation) << ",\n";

  return ss.str();
}

std::vector<cgtl::TaskAllocation> parse_task_allocations_mpsym(
  std::string const &task_allocations_str)
{ return split_task_allocations(task_allocations_str); }
