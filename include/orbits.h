#ifndef _GUARD_ORBITS_H
#define _GUARD_ORBITS_H

#include <algorithm>
#include <initializer_list>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "dump.h"
#include "perm_set.h"
#include "schreier_structure.h"

/**
 * @file orbits.h
 * @brief Free standing orbit calculation function(s).
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

class Orbit : public std::vector<unsigned>
{
  friend std::ostream &operator<<(std::ostream &os, Orbit const &o);

public:
  using const_iterator = std::vector<unsigned>::const_iterator;

  Orbit()
  : std::vector<unsigned>()
  {}

  Orbit(std::initializer_list<unsigned> orbit)
  : std::vector<unsigned>(orbit)
  {}

  template<typename IT>
  Orbit(IT first, IT last)
  : std::vector<unsigned>(first, last)
  {}

  static Orbit generate(unsigned x,
                        PermSet const &generators,
                        std::shared_ptr<SchreierStructure> ss = nullptr);

  bool generated_by(unsigned x,
                    PermSet const &generators) const;

  void update(PermSet const &generators_old,
              Perm const &generator_new,
              std::shared_ptr<SchreierStructure> ss = nullptr);

  bool contains(unsigned x) const
  { return std::find(_orbit.begin(), _orbit.end(), x) != _orbit.end(); }

private:
  void extend(PermSet const &generators,
              std::vector<unsigned> stack,
              std::unordered_set<unsigned> done,
              std::shared_ptr<SchreierStructure> ss);

  std::vector<unsigned> _orbit;
};

inline std::ostream &operator<<(std::ostream &os, Orbit const &o)
{
  os << DUMP(static_cast<std::vector<unsigned> const &>(o), "{}");
  return os;
}

class OrbitPartition
{
  friend std::ostream &operator<<(std::ostream &os, OrbitPartition const &op);

public:
  using const_iterator = std::vector<Orbit>::const_iterator;

  OrbitPartition(PermSet const &generators);

  unsigned num_partitions() const
  { return static_cast<unsigned>(_partitions.size()); }

  unsigned partition_index(unsigned x) const
  { return _partition_indices[x - 1u]; }

  Orbit const& operator[](unsigned i) const
  { return _partitions[i]; }

  const_iterator begin() const
  { return _partitions.begin(); }

  const_iterator end() const
  { return _partitions.end(); }

private:
  std::vector<Orbit> _partitions;
  std::vector<unsigned> _partition_indices;
};

inline std::ostream &operator<<(std::ostream &os, OrbitPartition const &op)
{
  os << DUMP(op._partitions, "{}");
  return os;
}

} // namespace cgtl

#endif // _GUARD_ORBITS_H
