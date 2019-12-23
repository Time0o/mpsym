#ifndef _GUARD_PROFILE_GENERATE_H
#define _GUARD_PROFILE_GENERATE_H

#include <algorithm>
#include <sstream>
#include <random>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

inline std::string generate_task_allocations(unsigned num_pes,
                                             unsigned num_tasks,
                                             unsigned num_task_allocations)
{
  throw std::logic_error("TODO: incomplete");

  static std::random_device rd{};
  static std::mt19937 m{rd()};

  std::stringstream ss;

  std::vector<unsigned> pes(num_pes);
  std::iota(pes.begin(), pes.end(), 1u);

  std::vector<unsigned> sample(num_tasks);
  for (unsigned i = 0u; i < num_task_allocations; ++i) {
    std::sample(pes.begin(), pes.end(), sample.begin(), num_tasks, m);

    ss << sample[0];
    for (unsigned j = 1u; j < num_tasks; ++j)
      ss << " " << sample[j];
    ss << "\n";
  }

  return ss.str();
}

#endif // _GUARD_PROFILE_GENERATE_H
