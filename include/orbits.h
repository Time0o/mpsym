#ifndef _GUARD_ORBITS_H
#define _GUARD_ORBITS_H

#include <vector>

#include "perm.h"

/**
 * @file orbits.h
 * @brief Free standing orbit calculation function(s).
 *
 * @author Timo Nicolai
 */

namespace cgtl
{

std::vector<unsigned>
orbit_of(unsigned x, std::vector<Perm> const &generators);

std::vector<std::vector<unsigned>>
orbit_partition(std::vector<Perm> const &generators);

} // namespace cgtl

#endif // _GUARD_ORBITS_H
