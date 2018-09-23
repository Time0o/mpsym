#ifndef _GUARD_SCHREIER_STRUCTURE_H
#define _GUARD_SCHREIER_STRUCTURE_H

#include <map>
#include <vector>

#include "perm.h"

/**
 * @file schreier_sims.h
 * @brief Schreier structure variant definitions.
 *
 * For details on the implementation difference between Schreier trees and
 * shallow Schreier trees consult \cite holt05.
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

struct SchreierStructure
{
  virtual void create_root(unsigned root) = 0;
  virtual void create_labels(std::vector<Perm> const &labels) = 0;
  virtual void create_edge(
    unsigned origin, unsigned destination, unsigned label) = 0;

  virtual unsigned root() const = 0;
  virtual std::vector<unsigned> nodes() const = 0;
  virtual std::vector<Perm> labels() const = 0;

  virtual bool contains(unsigned node) const = 0;
  virtual Perm transversal(unsigned origin) const = 0;
};

struct SchreierTree : public SchreierStructure
{
  SchreierTree(unsigned degree) : _degree(degree) {}

  void create_root(unsigned root) override { _root = root; }

  void create_labels(
    std::vector<Perm> const &labels) override { _labels = labels; }

  void create_edge(
    unsigned origin, unsigned destination, unsigned label) override {

    _edges[origin] = destination;
    _edge_labels[origin] = label;
  }

  unsigned root() const override { return _root; }

  std::vector<unsigned> nodes() const override {
    std::vector<unsigned> result {_root};

    for (auto const &node : _edges)
      result.push_back(node.first);

    return result;
  }

  std::vector<Perm> labels() const override { return _labels; }

  bool contains(unsigned node) const override {
    return (node == _root) || (_edges.find(node) != _edges.end());
  }

  Perm transversal(unsigned origin) const override
  {
    Perm result(_degree);

    unsigned current = origin;
    while(current != _root) {
      Perm const &label = _labels[_edge_labels.find(current)->second];
      result = label * result;
      current = _edges.find(current)->second;
    }

    return result;
  }

private:
  unsigned _degree;
  unsigned _root = 0;
  std::map<unsigned, unsigned> _edges;
  std::vector<Perm> _labels;
  std::map<unsigned, unsigned> _edge_labels;
};

} // namespace cgtl

#endif // _GUARD_SCHREIER_SIMS_H
