#include <chrono>
#include <cstdio>
#include <ctime>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "arch_graph.h"
#include "profile_utility.h"

#define RUNS 100

using cgtl::ArchGraph;
using cgtl::TaskMapping;

static long run(ArchGraph const &ag, unsigned num_tasks,
  ArchGraph::MappingVariant mapping_variant)
{
  static std::default_random_engine re(time(nullptr));

  std::uniform_int_distribution<unsigned> dist(0u, ag.num_processors() - 1u);

  std::vector<unsigned> task_mapping1(num_tasks);
  std::vector<unsigned> task_mapping2(num_tasks);

  long exec_ticks = 0;
  int exec_errors = 0;

  for (int run = 0; run < RUNS; ++run) {
    for (unsigned i = 0u; i < num_tasks; ++i) {
      task_mapping1[i] = dist(re);
      task_mapping2[i] = dist(re);
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    TaskMapping tm1 = ag.mapping(task_mapping1, mapping_variant);
    TaskMapping tm2 = ag.mapping(task_mapping2, mapping_variant);

    exec_ticks += std::chrono::duration_cast<std::chrono::microseconds>(
      std::chrono::high_resolution_clock::now() - t0).count();

    if (mapping_variant == ArchGraph::MAP_APPROX) {
      auto tm1_correct = ag.mapping(task_mapping1, ArchGraph::MAP_BRUTEFORCE);
      auto tm2_correct = ag.mapping(task_mapping2, ArchGraph::MAP_BRUTEFORCE);

      if (tm1.equivalence_class() != tm1_correct.equivalence_class())
        ++exec_errors;
      if (tm2.equivalence_class() != tm2_correct.equivalence_class())
        ++exec_errors;
    }
  }

  if (exec_errors > 0) {
    printf("WARNING: %d/%d approximate mappings incorrect\n",
           exec_errors, 2 * RUNS);
  }

  return exec_ticks;
}

int main()
{
/*
  std::vector<std::string> arch_graphs {
    resource_path("mcsoc.lua")
  };

  std::vector<std::pair<char const *, ArchGraph::MappingVariant>> variants {
    {"bruteforce", ArchGraph::MAP_BRUTEFORCE},
    {"approximation", ArchGraph::MAP_APPROX}
  };

  for (auto const &var : variants) {
    long exec_ticks = 0;

    for (auto const &path : arch_graphs) {
      ArchGraph ag;

      ag.fromlua(path);
      ag.complete();

      unsigned num_tasks = ag.num_processors();
      do {
        printf("Comparing mappings (%u task(s) on %u processors, %d mappings)\n",
               num_tasks, ag.num_processors(), RUNS);

        long ticks = run(ag, num_tasks, std::get<1>(var));
        printf("Execution time: %ld microseconds\n", ticks);

        exec_ticks += run(ag, num_tasks, std::get<1>(var));
      } while (num_tasks >>= 1u);
    }

    printf("Total execution time (%s): %ld microseconds\n",
           std::get<0>(var), exec_ticks);
  }
*/
}
