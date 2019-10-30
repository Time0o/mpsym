#ifndef _GUARD_SOLVABLE_BSGS_H
#define _GUARD_SOLVABLE_BSGS_H

#include <stdexcept>

#include "bsgs.h"
#include "perm_set.h"

namespace cgtl
{

typedef std::runtime_error BSGS_solve_error;

void solve_bsgs(BSGS &bsgs, PermSet const &generators);

} // namespace cgtl

#endif // _GUARD_SOLVABLE_BSGS_H
