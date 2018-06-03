#ifndef _GUARD_PERM_H
#define _GUARD_PERM_H

#include <map>
#include <vector>

namespace cgtl
{
class Perm
{
public:
  Perm(unsigned degree = 1);
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

  std::vector<unsigned> orbit() const;
  Perm transversal(unsigned origin) const;
  std::vector<Perm> transversals(std::vector<unsigned> const &origins) const;

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
  PermGroup(unsigned degree, std::vector<Perm> const &generators);

  static PermGroup symmetric(unsigned degree);
  static PermGroup cyclic(unsigned degree);
  static PermGroup alternating(unsigned degree);

  static std::vector<unsigned> orbit(unsigned alpha,
    std::vector<Perm> const &generators, SchreierTree &st);

  static std::pair<Perm, unsigned> strip(Perm const &perm,
    std::vector<unsigned> const &base, std::vector<Perm> const &generators,
    std::vector<SchreierTree> const &sts);

  static void schreier_sims(std::vector<unsigned> &base,
    std::vector<Perm> &generators, std::vector<SchreierTree> &schreier_trees);

  unsigned degree() const { return _n; }
  unsigned order() const { return _order; }
  std::vector<unsigned> base() const { return _base; };
  std::vector<Perm> strong_generating_set() const { return _strong_generating_set; };
  std::vector<SchreierTree> schreier_trees() const { return _schreier_trees; };

  bool is_member(Perm const &perm);

private:
  unsigned _n;
  unsigned _order;
  std::vector<unsigned> _base;
  std::vector<Perm> _strong_generating_set;
  std::vector<SchreierTree> _schreier_trees;
};

} // namespace cgtl

#endif // _GUARD_PERM_H
