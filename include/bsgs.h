#ifndef _GUARD_BSGS_H
#define _GUARD_BSGS_H

#include <utility>
#include <vector>

#include "perm.h"
#include "schreier_tree.h"

namespace cgtl
{

class BSGS
{
public:
  struct BaseElem {
    BaseElem(SchreierTree const &st) : _st(st) {}

    unsigned elem() const { return _st.root(); }
    std::vector<unsigned> orbit() const { return _st.nodes(); }
    Perm transversal(unsigned o) const { return _st.transversal(o); }
    std::vector<Perm> transversals() const;

  private:
    SchreierTree const &_st;
  };

  typedef std::vector<BaseElem>::size_type size_type;
  typedef std::vector<BaseElem>::const_iterator const_iterator;

  BSGS(std::vector<unsigned> const &base, std::vector<Perm> const &generators);
  BSGS(std::vector<Perm> const &generators) : BSGS({}, generators) {};

  const_iterator begin() const { return _base_elems.begin(); }
  const_iterator end() const { return _base_elems.end(); }
  BaseElem const& operator[](unsigned const i) const { return _base_elems[i]; }

  std::vector<unsigned> base() const { return _base; }
  std::vector<Perm> sgs() const { return _sgs; }
  std::vector<std::vector<unsigned>> orbits() const;
  size_type size() const { return _base_elems.size(); }
  std::pair<Perm, unsigned> strip(Perm const &perm) const;

private:
  std::vector<unsigned> _base;
  std::vector<Perm> _sgs;

  std::vector<SchreierTree> _schreier_trees;
  std::vector<BaseElem> _base_elems;
};

} // namespace cgtl

#endif // _GUARD_BSGS_H
