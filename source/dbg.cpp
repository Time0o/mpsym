#include "dbg.hpp"

#ifndef NDEBUG

namespace mpsym
{

namespace internal
{

namespace dbg
{

int Dbg::loglevel = WARN;
std::ostream *Dbg::out = &std::cout;

} // namespace dbg

} // namespace internal

} // namespace mpsym

#endif
