#ifndef _GUARD_PROFILE_ARGS_H
#define _GUARD_PROFILE_ARGS_H

#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <vector>

class VariantOption
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
  std::vector<char const *>::size_type choice_index(char const *choice) const
  {
    for (auto i = 0u; i < _choices.size(); ++i) {
      if (strcmp(_choices[i], choice) == 0)
        return i;
    }

    throw std::invalid_argument("invalid parameter choice");
  }

  std::vector<char const *> _choices;
  std::vector<char const *>::size_type _current_choice;
};

#endif // _GUARD_PROFILE_ARGS_H
