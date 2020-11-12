#ifndef GUARD_ORBITS_H
#define GUARD_ORBITS_H

#include <algorithm>
#include <initializer_list>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "dump.hpp"

namespace mpsym
{

namespace internal
{

class PermSet;
class SchreierStructure;

class Orbit
{
  friend std::ostream &operator<<(std::ostream &os, Orbit const &o);

public:
  using value_type = unsigned;
  using size_type = std::vector<unsigned>::size_type;

  using const_iterator = std::vector<unsigned>::const_iterator;

  Orbit()
  {}

  Orbit(std::initializer_list<unsigned> elements)
  : _elements(elements)
  {}

  template<typename IT>
  Orbit(IT first, IT last)
  : _elements(first, last)
  {}

  static Orbit generate(unsigned x,
                        PermSet const &generators,
                        std::shared_ptr<SchreierStructure> ss = nullptr);

  bool operator==(Orbit const &other) const;

  bool operator!=(Orbit const &other) const
  { return !(*this == other); }

  bool generated_by(unsigned x,
                    PermSet const &generators) const;

  void update(PermSet const &generators_old,
              PermSet const &generators_new,
              std::shared_ptr<SchreierStructure> ss = nullptr);

  void insert(unsigned x)
  { _elements.push_back(x); }

  template<typename IT>
  bool erase(unsigned x)
  {
    for (auto it = begin(); it != end(); ++it) {
      if (*it == x) {
        erase(it);
        return true;
      }
    }

    return false;
  }

  template<typename IT>
  IT erase(IT it)
  { return _elements.erase(it); }

  bool empty() const
  { return _elements.empty(); }

  size_type size() const
  { return _elements.size(); }

  const_iterator begin() const
  { return _elements.begin(); }

  const_iterator end() const
  { return _elements.end(); }

  bool contains(unsigned x) const
  { return std::find(begin(), end(), x) != end(); }

private:
  void extend(PermSet const &generators,
              std::vector<unsigned> stack,
              std::unordered_set<unsigned> done,
              std::shared_ptr<SchreierStructure> ss);

  std::vector<unsigned> _elements;
};

inline std::ostream &operator<<(std::ostream &os, Orbit const &orbit)
{
  os << DUMP_CUSTOM(orbit._elements);
  return os;
}

class OrbitPartition
{
  friend std::ostream &operator<<(std::ostream &os, OrbitPartition const &op);

public:
  using const_iterator = std::vector<Orbit>::const_iterator;

  explicit OrbitPartition(unsigned degree);

  OrbitPartition(unsigned degree, std::vector<Orbit> const &partitions);

  OrbitPartition(unsigned degree, std::vector<int> const &partition_indices);

  OrbitPartition(unsigned degree, PermSet const &generators);

  bool operator==(OrbitPartition const &other) const
  { return _partition_indices == other._partition_indices; }

  bool operator!=(OrbitPartition const &other) const
  { return !(*this == other); }

  std::vector<OrbitPartition> split(OrbitPartition const &split) const;

  unsigned num_partitions() const
  { return static_cast<unsigned>(_partitions.size()); }

  int partition_index(unsigned x) const
  { return _partition_indices[x - 1u]; }

  void remove_from_partition(unsigned x);

  void change_partition(unsigned x, int i);

  Orbit const& operator[](unsigned i) const
  { return _partitions[i]; }

  const_iterator begin() const
  { return _partitions.begin(); }

  const_iterator end() const
  { return _partitions.end(); }

private:
  void add_to_partition(unsigned x, int i);
  void update_partitions();
  void update_partition_indices();

  std::vector<Orbit> _partitions;
  std::vector<int> _partition_indices;
};

inline std::ostream &operator<<(std::ostream &os, OrbitPartition const &op)
{
  os << DUMP_CUSTOM(op._partitions, "{}");
  return os;
}

} // namespace internal

} // namespace mpsym

#endif // GUARD_ORBITS_H
