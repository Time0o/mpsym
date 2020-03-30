#ifndef _GUARD_ARCH_GRAPH_SYSTEM_H
#define _GUARD_ARCH_GRAPH_SYSTEM_H

#include <memory>
#include <string>

#include "bsgs.h"
#include "partial_perm_inverse_semigroup.h"
#include "perm_group.h"
#include "task_mapping.h"
#include "task_orbits.h"

namespace mpsym
{

class ArchGraphSystem
{
public:
  using AutomorphismOptions = BSGS::Options;

  enum class ReprMethod {
    ITERATE,
    LOCAL_SEARCH,
    ORBITS,
    AUTO = ITERATE
  };

  struct ReprOptions {
    static ReprOptions fill_defaults(ReprOptions const *options)
    {
      static ReprOptions default_options;
      return options ? *options : default_options;
    }

    ReprMethod method = ReprMethod::AUTO;
    unsigned offset = 0u;
    bool match = true;
  };

  static std::shared_ptr<ArchGraphSystem> from_lua(std::string const &lua);

  virtual std::string to_gap() const = 0;

  virtual unsigned num_processors() const
  { throw std::logic_error("not implemented"); }

  virtual unsigned num_channels() const
  { throw std::logic_error("not implemented"); }

  bool automorphisms_ready() const
  { return _automorphisms_valid; }

  void reset_automorphisms()
  { _automorphisms_valid = false; }

  BSGS::order_type num_automorphisms(AutomorphismOptions const *options = nullptr)
  { return num_automorphisms_(options); }

  PermGroup automorphisms(AutomorphismOptions const *options = nullptr)
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

private:
  virtual BSGS::order_type num_automorphisms_(AutomorphismOptions const *options)
  { return automorphisms(options).order(); }

  virtual PermGroup automorphisms_(AutomorphismOptions const *options) = 0;

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
                                    TaskOrbits *orbits,
                                    ReprOptions const *options);

  PermGroup _automorphisms;
  bool _automorphisms_valid = false;
};

// TODO: outstream operator

} // namespace mpsym

#endif // _GUARD_ARCH_GRAPH_SYSTEM_H
