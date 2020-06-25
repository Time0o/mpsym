#ifndef _GUARD_PERMLIB_H
#define _GUARD_PERMLIB_H

#include <iterator>

namespace boost
{

// permlib needs this to compile and I couldn't be bothered to build this
// program against an appropriate version of boost instead
template<typename ForwardIt>
ForwardIt next(
  ForwardIt it,
  typename std::iterator_traits<ForwardIt>::difference_type n = 1) {

  return std::next(it, n);
}

}

#include <permlib/permlib_api.h>

#include <permlib/construct/random_schreier_sims_construction.h>
#include <permlib/generator/bsgs_random_generator.h>
#include <permlib/permutation.h>
#include <permlib/transversal/explicit_transversal.h>
#include <permlib/transversal/shallow_schreier_tree_transversal.h>

#endif // _GUARD_PERMLIB_H
