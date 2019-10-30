#ifndef _GUARD_ORBITS_H
#define _GUARD_ORBITS_H

#include <vector>

#include "perm_set.h"

/**
 * @file orbits.h
 * @brief Free standing orbit calculation function(s).
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

std::vector<unsigned>
orbit_of(unsigned x, PermSet const &generators);

std::vector<std::vector<unsigned>>
orbit_partition(PermSet const &generators);

} // namespace cgtl

#endif // _GUARD_ORBITS_H
