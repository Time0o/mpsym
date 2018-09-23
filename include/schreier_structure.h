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

/** Data structure for space efficient storage of orbit transversals.
 *
 * The *orbit* of a group element \f$g \in G\f$ where \f$G\f$ acts on the set
 * \f$\Omega\f$ is the set \f$G(x) = \{x^g \mid g \in G\}\f$. The *orbit
 * transversals* corresponding to this orbit are \f$\{u_{\beta} \mid \beta \in
 * G(x)\}\f$ where \f$u_{\beta} \in G\f$ maps \f$x\f$ to \f$\beta\f$. Instead
 * of storing these explicitly (which can might consume a lot of space for
 * large orbits (possible for large sets \f$\Omega\f$), this data structure
 * only stores a spanning tree for the orbit graph connecting every element
 * \f$\beta' \in G(x)\f$ to every other element \f$\beta\f$ via a directed edge
 * labelled with \f$g \in G(x)\f$ if \f$\exists g \in G \colon \beta^g = \beta'\f$.
 * The spanning tree is rooted at \f$x\f$.
 */
struct SchreierTree : public SchreierStructure
{
  SchreierTree(unsigned degree) : _degree(degree) {}

  /** Set a spanning tree's root.
   *
   * \param root
   *     the element \f$x\f$ as per the notation used in this struct's
   *     definition
   */
  void create_root(unsigned root) override { _root = root; }

  /** Set a Schreier tree's possible edge labels.
   *
   * This function must be called before the first call to edge().
   *
   * \param labels possible permutation edge labels
   */
  void create_labels(
    std::vector<Perm> const &labels) override { _labels = labels; }

  /** Add an edge to a Schreier tree.
   *
   * \param origin
   *     an element \f$\beta' \in \Omega\f$ as per the notation used in this
   *     struct's definition
   *
   * \param destination
   *     an element \f$\beta \in \Omega\f$ as per the notation used in this
   *     struct's definition
   *
   * \param perm
   *     a perm object describing a permutation \f$g \in G\f$ with \f$\beta^g =
   *     \beta'\f$
   */
  void create_edge(
    unsigned origin, unsigned destination, unsigned label) override {

    _edges[origin] = destination;
    _edge_labels[origin] = label;
  }

  /** Obtain a Schreier tree's root element.
   *
   * \return the Schreier tree's root, i.e. the element \f$x\f$ as per the
   *         notation used in this struct's definition
   */
  unsigned root() const override { return _root; }

  /** Obtain a Schreier tree's nodes.
   *
   * \return the Schreier tree's nodes, i.e. the orbit \f$G(x)\f$ as per the
   *         notation used in this struct's definition, in the form of a vector
   *         with no particular guaranteed order
   */
  std::vector<unsigned> nodes() const override {
    std::vector<unsigned> result {_root};

    for (auto const &node : _edges)
      result.push_back(node.first);

    return result;
  }

  /** Obtain all of a Schreier tree's possible edge labels.
   *
   * \return this Schreier tree's unique edge labels in not particular order,
   *         in the form of a vector
   */
  std::vector<Perm> labels() const override { return _labels; }

  /** Check whether an element of \f$\Omega\f$ is a node of a Schreier tree.
   *
   * This is equivalent to testing whether the element is contained in the orbit
   * \f$G(x)\f$.
   *
   * \param node node to be testing for membership in this Schreier tree
   *
   * \return `true` if `node` is a node of this Schreier tree, else `false`
   */
  bool contains(unsigned node) const override {
    return (node == _root) || (_edges.find(node) != _edges.end());
  }

  /** Obtain an orbit transversal from a Schreier tree.
   *
   * \param origin
   *     an element \f$\beta\f$ of \f$\Omega\f$, this must be a node of this
   *     Schreier tree, otherwise this function's behaviour is undefined
   *
   * \return a permutation \f$g \in G\f$ such that \f$x^g = \beta\f$
   */
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
