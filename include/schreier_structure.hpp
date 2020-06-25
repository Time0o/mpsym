#ifndef GUARD_SCHREIER_STRUCTURE_H
#define GUARD_SCHREIER_STRUCTURE_H

#include <ostream>
#include <vector>

#include "perm.hpp"
#include "perm_set.hpp"

namespace mpsym
{

namespace internal
{

struct SchreierStructure
{
  friend std::ostream &operator<<(std::ostream &os, SchreierStructure const &ss);

  virtual void add_label(Perm const &label) = 0;

  virtual void create_edge(unsigned origin,
                           unsigned destination,
                           unsigned label) = 0;

  virtual unsigned root() const = 0;
  virtual std::vector<unsigned> nodes() const = 0;
  virtual PermSet labels() const = 0;

  virtual bool contains(unsigned node) const = 0;
  virtual bool incoming(unsigned node, Perm const &edge) const = 0;
  virtual Perm transversal(unsigned origin) const = 0;

private:
  virtual void dump(std::ostream& os) const = 0;
};

inline std::ostream &operator<<(std::ostream &os, SchreierStructure const &ss)
{
  ss.dump(os);
  return os;
}

} // namespace internal

} // namespace mpsym

#endif // GUARD_SCHREIER_SIMS_H
