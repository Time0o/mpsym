#ifndef _GUARD_PROFILE_PARSE_H
#define _GUARD_PROFILE_PARSE_H

#include <string>
#include <tuple>
#include <vector>

#include "perm.h"
#include "permlib.h"

std::tuple<unsigned, unsigned, std::string> parse_group(
  std::string const &group_str);

std::string parse_generators_gap(
  std::string const &gen_str);

std::vector<cgtl::Perm> parse_generators_mpsym(
  std::string const &gen_str);

std::vector<permlib::Permutation::ptr> parse_generators_permlib(
  std::string const &gen_str);

std::string parse_task_allocations_gap(
  std::string const &task_allocations_str);

#endif // _GUARD_PROFILE_PARSE_H
