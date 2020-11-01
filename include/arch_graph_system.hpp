#ifndef GUARD_ARCH_GRAPH_SYSTEM_H
#define GUARD_ARCH_GRAPH_SYSTEM_H

#include <memory>
#include <string>
#include <tuple>
#include <unordered_set>

#include "bsgs.hpp"
#include "perm_group.hpp"
#include "string.hpp"
#include "task_mapping.hpp"
#include "task_orbits.hpp"
#include "timeout.hpp"

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

  static std::shared_ptr<ArchGraphSystem> from_lua(
    std::string const &lua,
    std::vector<std::string> const &args = {});

  static std::shared_ptr<ArchGraphSystem> from_lua_file(
    std::string const &lua_file,
    std::vector<std::string> const &args = {})
  { return from_lua(util::read_file(lua_file), args); }

  static std::shared_ptr<ArchGraphSystem> from_json(
    std::string const &json);

  static std::shared_ptr<ArchGraphSystem> from_json_file(
    std::string const &json_file)
  { return from_json(util::read_file(json_file)); }

  virtual std::string to_gap() const = 0;
  virtual std::string to_json() const = 0;

  std::shared_ptr<ArchGraphSystem> expand_automorphisms() const;

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
    AutomorphismOptions const *options = nullptr,
    internal::timeout::aborted_type aborted = internal::timeout::unaborted())
  { return num_automorphisms_(options, aborted); }

  internal::PermGroup automorphisms(
    AutomorphismOptions const *options = nullptr,
    internal::timeout::aborted_type aborted = internal::timeout::unaborted())
  {
    if (!automorphisms_ready()) {
      _automorphisms = automorphisms_(options, aborted);
      _automorphisms_valid = true;
    }

    return _automorphisms;
  }

  void init_repr(
    AutomorphismOptions const *options = nullptr,
    internal::timeout::aborted_type aborted = internal::timeout::unaborted())
  {
    if (!repr_ready_())
      init_repr_(options, aborted);
  }

  bool repr_ready() const
  { return repr_ready_(); }

  void reset_repr()
  { reset_repr_(); }

  TaskMapping repr(
    TaskMapping const &mapping,
    ReprOptions const *options = nullptr,
    internal::timeout::aborted_type aborted = internal::timeout::unaborted())
  {
    if (!repr_ready_())
      init_repr();

    return repr_(mapping, options, nullptr, aborted);
  }

  std::tuple<TaskMapping, bool, unsigned> repr(
    TaskMapping const &mapping,
    TaskOrbits &orbits,
    ReprOptions const *options = nullptr,
    internal::timeout::aborted_type aborted = internal::timeout::unaborted())
  {
    if (!repr_ready_())
      init_repr();

    auto representative(repr_(mapping, options, &orbits, aborted));

    auto ins(orbits.insert(representative));

    return std::make_tuple(representative, ins.first, ins.second);
  }

  std::vector<TaskMapping> orbit(TaskMapping const &mapping);

private:
  virtual internal::BSGS::order_type num_automorphisms_(
    AutomorphismOptions const *options,
    internal::timeout::aborted_type aborted)
  { return automorphisms(options, aborted).order(); }

  virtual internal::PermGroup automorphisms_(
    AutomorphismOptions const *options,
    internal::timeout::aborted_type aborted) = 0;

  virtual void init_repr_(AutomorphismOptions const *,
                          internal::timeout::aborted_type )
  {}

  virtual TaskMapping repr_(TaskMapping const &mapping,
                            ReprOptions const *options,
                            TaskOrbits *orbits,
                            internal::timeout::aborted_type aborted);

  virtual bool repr_ready_() const
  { return automorphisms_ready(); }

  virtual void reset_repr_()
  { reset_automorphisms(); }

  bool automorphisms_symmetric(ReprOptions const *options);

  TaskMapping repr_symmetric(TaskMapping const &mapping,
                             ReprOptions const *options);

  static bool is_repr(TaskMapping const &tasks,
                      ReprOptions const *options,
                      TaskOrbits *orbits)
  {
    if (!options->match || !orbits)
      return false;

    return orbits->is_repr(tasks);
  }

  TaskMapping min_elem_iterate(TaskMapping const &tasks,
                               ReprOptions const *options,
                               TaskOrbits *orbits,
                               internal::timeout::aborted_type aborted);

  TaskMapping min_elem_orbits(TaskMapping const &tasks,
                              ReprOptions const *options,
                              TaskOrbits *orbits,
                              internal::timeout::aborted_type aborted);

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

} // namespace mpsym

#endif // GUARD_ARCH_GRAPH_SYSTEM_H
