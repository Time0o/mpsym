#ifndef _GUARD_SCHREIER_STRUCTURE_H
#define _GUARD_SCHREIER_STRUCTURE_H

#include <map>
#include <vector>

#include "perm.h"
#include "perm_set.h"

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
  virtual void create_labels(PermSet const &labels) = 0;
  virtual void create_edge(
    unsigned origin, unsigned destination, unsigned label) = 0;

  virtual unsigned root() const = 0;
  virtual std::vector<unsigned> nodes() const = 0;
  virtual PermSet labels() const = 0;

  virtual bool contains(unsigned node) const = 0;
  virtual bool incoming(unsigned node, Perm const &edge) const = 0;
  virtual Perm transversal(unsigned origin) const = 0;
};

struct ExplicitTransversals : public SchreierStructure
{
  ExplicitTransversals(unsigned degree) : _degree(degree) {}

  void create_root(unsigned root) override;
  void create_labels(PermSet const &labels) override;
  void create_edge(
    unsigned origin, unsigned destination, unsigned label) override;

  unsigned root() const override;
  std::vector<unsigned> nodes() const override;
  PermSet labels() const override;

  bool contains(unsigned node) const override;
  bool incoming(unsigned node, Perm const &edge) const override;
  Perm transversal(unsigned origin) const override;

private:
  unsigned _degree;
  unsigned _root = 0;
  PermSet _labels;
  std::map<unsigned, Perm> _orbit;
};

struct SchreierTree : public SchreierStructure
{
  SchreierTree(unsigned degree) : _degree(degree) {}

  void create_root(unsigned root) override;
  void create_labels(PermSet const &labels) override;
  void create_edge(
    unsigned origin, unsigned destination, unsigned label) override;

  unsigned root() const override;
  std::vector<unsigned> nodes() const override;
  PermSet labels() const override;

  bool contains(unsigned node) const override;
  bool incoming(unsigned node, Perm const &edge) const override;
  Perm transversal(unsigned origin) const override;

private:
  unsigned _degree;
  unsigned _root = 0;
  std::map<unsigned, unsigned> _edges;
  PermSet _labels;
  std::map<unsigned, unsigned> _edge_labels;
};

} // namespace cgtl

#endif // _GUARD_SCHREIER_SIMS_H
