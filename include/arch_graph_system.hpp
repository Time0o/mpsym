#ifndef GUARD_ARCH_GRAPH_SYSTEM_H
#define GUARD_ARCH_GRAPH_SYSTEM_H

#include <memory>
#include <string>
#include <unordered_set>

#include "bsgs.hpp"
#include "perm_group.hpp"
#include "task_mapping.hpp"
#include "task_orbits.hpp"

namespace mpsym
{

using AutomorphismOptions = internal::BSGSOptions;

struct ReprOptions
{
  enum class Method {
    ITERATE,
    LOCAL_SEARCH,
    ORBITS,
    AUTO = ITERATE
  };

  enum class Variant {
    LOCAL_SEARCH_BFS,
    LOCAL_SEARCH_DFS,
    LOCAL_SEARCH_SA_LINEAR
  };

  static ReprOptions fill_defaults(ReprOptions const *options)
  {
    static ReprOptions default_options;
    return options ? *options : default_options;
  }

  Method method = Method::AUTO;
  Variant variant = Variant::LOCAL_SEARCH_BFS;

  unsigned offset = 0u;

  bool match = true;
  bool optimize_symmetric = true;

  bool local_search_invert_generators = false;
  unsigned local_search_append_generators = 0u;
  unsigned local_search_sa_iterations = 100u;
  double local_search_sa_T_init = 1.0;
};

class ArchGraphSystem
{
public:
  virtual ~ArchGraphSystem() = default;

  static std::shared_ptr<ArchGraphSystem> from_lua_file(
    std::string const &lua_file,
    std::vector<std::string> const &args = {});

  static std::shared_ptr<ArchGraphSystem> from_lua(
    std::string const &lua,
    std::vector<std::string> const &args = {});

  virtual std::string to_gap() const = 0;

  virtual std::string to_json();

  virtual unsigned num_processors() const
  { throw std::logic_error("not implemented"); }

  virtual unsigned num_channels() const
  { throw std::logic_error("not implemented"); }

  bool automorphisms_ready() const
  { return _automorphisms_valid; }

  void reset_automorphisms()
  {
    _automorphisms_valid = false;
    _automorphisms_is_shifted_symmetric_valid = false;
  }

  internal::BSGS::order_type num_automorphisms(
    AutomorphismOptions const *options = nullptr)
  { return num_automorphisms_(options); }

  unsigned num_automorphism_orbits(
    unsigned num_tasks,
    bool unique_tasks,
    AutomorphismOptions const *options = nullptr);

  std::vector<unsigned> automorphism_orbit_sizes(
    unsigned num_tasks,
    bool unique_tasks,
    AutomorphismOptions const *options = nullptr);

  internal::PermGroup automorphisms(
    AutomorphismOptions const *options = nullptr)
  {
    if (!automorphisms_ready()) {
      _automorphisms = automorphisms_(options);
      _automorphisms_valid = true;
    }

    return _automorphisms;
  }

  void init_repr(AutomorphismOptions const *options = nullptr)
  {
    if (!repr_ready_())
      init_repr_(options);
  }

  bool repr_ready() const
  { return repr_ready_(); }

  void reset_repr()
  { reset_repr_(); }

  TaskMapping repr(TaskMapping const &mapping,
                   TaskOrbits *orbits = nullptr,
                   ReprOptions const *options = nullptr)
  {
    if (!repr_ready_())
      init_repr();

    return repr_(mapping, orbits, options);
  }

  std::vector<TaskMapping> orbit(TaskMapping const &mapping);

private:
  virtual internal::BSGS::order_type num_automorphisms_(
    AutomorphismOptions const *options)
  { return automorphisms(options).order(); }

  virtual internal::PermGroup automorphisms_(
    AutomorphismOptions const *options) = 0;

  virtual void init_repr_(AutomorphismOptions const *options)
  { automorphisms(options); }

  virtual bool repr_ready_() const
  { return automorphisms_ready(); }

  virtual void reset_repr_()
  { reset_automorphisms(); }

  virtual TaskMapping repr_(TaskMapping const &mapping,
                            TaskOrbits *orbits,
                            ReprOptions const *options);

  static bool is_repr(TaskMapping const &tasks,
                      TaskOrbits *orbits,
                      ReprOptions const *options)
  {
    if (!options->match || !orbits)
      return false;

    return orbits->is_repr(tasks);
  }

  TaskMapping min_elem_iterate(TaskMapping const &tasks,
                               TaskOrbits *orbits,
                               ReprOptions const *options);

  TaskMapping min_elem_orbits(TaskMapping const &tasks,
                              TaskOrbits *orbits,
                              ReprOptions const *options);

  TaskMapping min_elem_local_search(TaskMapping const &tasks,
                                    ReprOptions const *options);

  static internal::PermSet local_search_augment_gens(
    internal::PermGroup const &automorphisms,
    ReprOptions const *options);

  TaskMapping min_elem_local_search_sa(TaskMapping const &tasks,
                                       unsigned task_min,
                                       unsigned task_max,
                                       ReprOptions const *options);

  static double local_search_sa_schedule_T(unsigned i,
                                           ReprOptions const *options);

  static double local_search_sa_value(TaskMapping const &representative,
                                      unsigned task_min,
                                      unsigned task_max);

  TaskMapping min_elem_symmetric(TaskMapping const &tasks,
                                 unsigned task_min,
                                 unsigned task_max,
                                 ReprOptions const *options);

  TaskMapping random_task_mapping(unsigned num_tasks, bool unique_tasks);
  std::unordered_set<TaskMapping> task_mapping_orbit(TaskMapping const &tasks);

  internal::PermGroup _automorphisms;
  bool _automorphisms_valid = false;

  bool _automorphisms_is_shifted_symmetric;
  bool _automorphisms_is_shifted_symmetric_valid = false;

  unsigned _automorphisms_smp;
  unsigned _automorphisms_lmp;
};

// TODO: outstream operator

} // namespace mpsym

#endif // GUARD_ARCH_GRAPH_SYSTEM_H
