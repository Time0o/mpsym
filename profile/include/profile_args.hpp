#ifndef _GUARD_PROFILE_ARGS_H
#define _GUARD_PROFILE_ARGS_H

#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <unordered_set>
#include <vector>
#include <iostream>

#include "util.hpp"

namespace profile
{

struct VariantOptionBase
{
  std::vector<char const *>::size_type choice_index(char const *choice) const
  {
    for (auto i = 0u; i < _choices.size(); ++i) {
      if (strcmp(_choices[i], choice) == 0)
        return i;
    }

    std::string msg("invalid parameter choice: ");

    throw std::invalid_argument(msg + choice);
  }

  std::vector<char const *> _choices;
};

class VariantOption : private VariantOptionBase
{
public:
  VariantOption(std::initializer_list<char const *> choices)
  {
    _choices.emplace_back("unset");
    _choices.insert(_choices.end(), choices.begin(), choices.end());

    _current_choice = 0u;
  }

  void set(char const *choice)
  { _current_choice = choice_index(choice); }

  char const *get() const
  { return _choices[_current_choice]; }

  bool is_set() const
  { return _current_choice != 0u; }

  bool is(char const *choice) const
  { return _current_choice == choice_index(choice); }

private:
  std::vector<char const *>::size_type _current_choice;
};

class VariantOptionSet : private VariantOptionBase
{
public:
  VariantOptionSet(std::initializer_list<char const *> choices)
  {
    _choices.emplace_back("unset");
    _choices.insert(_choices.end(), choices.begin(), choices.end());
  }

  void set(char const *choice)
  {
    if (!is_set(choice))
      _selected.insert(choice_index(choice));
  }

  void unset(char const *choice)
  { _selected.erase(choice_index(choice)); }

  std::vector<char const *> get() const
  {
    std::vector<char const *> res;

    for (auto i : _selected)
      res.push_back(_choices[i]);

    return res;
  }

  bool is_set(char const *choice) const
  { return _selected.find(choice_index(choice)) != _selected.end(); }

private:
  std::unordered_set<std::vector<char const *>::size_type> _selected;
};

template<typename FUNC>
void foreach_option(char const *optarg, FUNC &&f)
{
  for (auto const &option : mpsym::util::split(optarg, ",")) {
    if (!option.empty())
      f(option);
  }
}

#define CHECK_OPTION(cond, msg) \
  if (!(cond)) { \
    usage(std::cerr); \
    error(msg); \
    return EXIT_FAILURE; \
  }

#define CHECK_ARGUMENT(arg) \
  if (optind == argc) { \
    usage(std::cerr); \
    error(arg, "argument is mandatory"); \
    return EXIT_FAILURE; \
  }

#define OPEN_STREAM(var, arg) \
  try { \
    var.open(arg); \
  } catch (std::runtime_error const &e) { \
    error(e.what(), ":", arg); \
    return EXIT_FAILURE; \
  }

} // namespace profile

#endif // _GUARD_PROFILE_ARGS_H
