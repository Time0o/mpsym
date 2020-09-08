#ifndef GUARD_RANDOM_H
#define GUARD_RANDOM_H

#include <random>

namespace mpsym
{

namespace util
{

inline std::mt19937 random_engine()
{ return std::mt19937{std::random_device{}()}; }

} // namespace util

} // namespace mpsym

#endif // GUARD_RANDOM_H
