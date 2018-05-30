#ifndef _GUARD_PERM_H
#define _GUARD_PERM_H

#include <cassert>
#include <initializer_list>
#include <vector>

#ifndef NDEBUG
#include <set>
#endif

namespace CGTL
{

template<unsigned N> class Perm;

template <unsigned M1, unsigned M2, unsigned MAX = std::max(M1, M2)>
Perm<MAX> operator*(Perm<M1> const &lhs, Perm<M2> const &rhs) {
  Perm<MAX> result;

  if (MAX == M1) {
    for (unsigned i = 0u; i < M2; ++i)
      result._perm[i] = lhs[rhs[i + 1]];

    for(unsigned j = M2; j < M1; ++j)
       result._perm[j] = lhs[j + 1];
  } else {
    for (unsigned i = 0u; i < M2; ++i) {
      unsigned tmp = rhs[i + 1];
      result._perm[i] = (tmp <= M1) ? lhs[tmp] : tmp;
    }
  }

  return result;
}

template<unsigned N>
class Perm
{
template <unsigned M1, unsigned M2, unsigned MAX>
friend Perm<MAX> operator*(Perm<M1> const &lhs, Perm<M2> const &rhs);

public:
  Perm() : _perm(N) {
    for (unsigned i = 0u; i < N; ++i)
       _perm[i] = i + 1;
  }

  template<unsigned M>
  Perm(Perm<M> const &old) {
    static_assert(N >= M, "can not narrow down permutation");

    for (unsigned i = 1u; i <= std::min(N, M); ++i)
      _perm.push_back(old[i]);

    if (N > M) {
      for (unsigned j = M + 1; j <= N; ++j)
        _perm.push_back(j);
    }
  }

  Perm(std::vector<unsigned> cycle) : Perm() {
#ifndef NDEBUG
    assert(("cycle implausibly long", cycle.size() <= N));

    {
      std::set<unsigned> tmp(cycle.begin(), cycle.end());

      assert(("cycle contains element > N", *tmp.rbegin() <= N));
      assert(("cycle contains duplicate elements", tmp.size() == cycle.size()));
    }
#endif
    for (size_t i = 1u; i < cycle.size(); ++i) {
       (*this)[cycle[i - 1]] = cycle[i];
    }
    (*this)[cycle.back()] = cycle[0];
  }

  Perm(std::vector<std::vector<unsigned>> cycles) {
    Perm<N> result(cycles.back());
    for (auto i = cycles.rbegin() + 1; i != cycles.rend(); ++i)
      result = Perm(*i) * result;

    _perm = result._perm;
  }

  unsigned const& operator[](unsigned const i) const {
    assert(("permutation index out of range", i > 0u && i <= N));
    return _perm[i - 1];
  }

  unsigned& operator[](unsigned const i) {
    return const_cast<unsigned&>(static_cast<const Perm*>(this)->operator[](i));
  }

private:
  std::vector<unsigned> _perm;
};

}

#endif // _GUARD_PERM_H
