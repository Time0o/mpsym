#ifndef _GUARD_PERM_H
#define _GUARD_PERM_H

#include <map>
#include <vector>

namespace cgtl
{

class Perm
{
public:
  Perm() : _n(0), _perm(0) {};
  Perm(unsigned degree);
  Perm(std::vector<unsigned> const &perm);
  Perm(unsigned n, std::vector<std::vector<unsigned>> const &cycles);

  unsigned const& operator[](unsigned const i) const;
  Perm operator~() const;
  bool operator==(Perm const &rhs) const;
  bool operator!=(Perm const &rhs) const;
  Perm& operator*=(Perm const &rhs);

  unsigned degree() const { return _n; }
  bool id() const;

private:
  unsigned _n;
  std::vector<unsigned> _perm;
};

std::ostream& operator<<(std::ostream& stream, const Perm &perm);
Perm operator*(Perm const &lhs, Perm const &rhs);

struct SchreierTree
{
  SchreierTree(unsigned degree) : _degree(degree) {}

  void create_root(unsigned root) { _root = root; }

  void create_edge(unsigned origin, unsigned destination, Perm const &perm) {
    _edges[origin] = destination;
    _labels[origin] = perm;
  }

  Perm transversal(unsigned origin) const;

  bool contains(unsigned node) const {
    return (node == _root) || (_edges.find(node) != _edges.end());
  }

private:
  unsigned _degree;
  unsigned _root;
  std::map<unsigned, unsigned> _edges;
  std::map<unsigned, Perm> _labels;
};

class PermGroup
{
public:
  PermGroup(unsigned degree, std::vector<Perm> const &generators)
    : _n(degree), _generators(generators) {}

  static std::vector<unsigned> orbit(unsigned alpha,
    std::vector<Perm> const &generators, SchreierTree &st);

  static std::pair<Perm, unsigned> strip(Perm const &perm,
    std::vector<unsigned> const &base, std::vector<Perm> const &generators,
    std::vector<SchreierTree> const &sts);

  static void schreier_sims(std::vector<unsigned> &base,
    std::vector<Perm> &generators);

  unsigned degree() const { return _n; }

private:
  unsigned _n;
  std::vector<Perm> _generators;
};

} // namespace cgtl

#endif // _GUARD_PERM_H
