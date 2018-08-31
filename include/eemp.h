#ifndef _GUARD_EEMP_H
#define _GUARD_EEMP_H

#include <vector>

#include "partial_perm.h"

namespace cgtl
{

class EEMP
{
public:
  static void action_components(
    std::vector<unsigned> const &alpha,
    std::vector<PartialPerm> const &generators);
};

} // namespace cgtl

#endif // _GUARD_EEMP_H
