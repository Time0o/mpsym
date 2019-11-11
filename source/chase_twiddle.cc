/* This implementation of Chase's Twiddle is based on one by Matthew Belmonte,
 * whose original copyright licence is included below. Changes to the
 * implementation are purely cosmetical in order to comply with the coding
 * standard of this project.
 */
/*
 * Coded by Matthew Belmonte <mkb4@Cornell.edu>, 23 March 1996.  This
 * implementation Copyright (c) 1996 by Matthew Belmonte.  Permission for use and
 * distribution is hereby granted, subject to the restrictions that this
 * copyright notice and reference list be included in its entirety, and that any
 * and all changes made to the program be clearly noted in the program text.
 *
 * This software is provided 'as is', with no warranty, express or implied,
 * including but not limited to warranties of merchantability or fitness for a
 * particular purpose.  The user of this software assumes liability for any and
 * all damages, whether direct or consequential, arising from its use.  The
 * author of this implementation will not be liable for any such damages.
 *
 * Reference:
 *
 * Phillip J Chase, `Algorithm 382: Combinations of M out of N Objects [G6]',
 * Communications of the Association for Computing Machinery 13:6:368 (1970).
 *
 * The returned indices x, y, and z in this implementation are decremented by 1,
 * in order to conform to the C language array reference convention.  Also, the
 * parameter 'done' has been replaced with a Boolean return value.
 */
#include <vector>

#include "chase_twiddle.h"

namespace cgtl
{

void ChaseTwiddle::init(int n, int k, std::vector<int> &p)
{
  int i;

  p[0] = n + 1;
  for (i = 1u; i <= n - k; ++i)
    p[i] = 0;

  while (i <= n) {
    p[i] = i + k - n;
    ++i;
  }

  p[n + 1] = -2;

  if (k == 0)
    p[1] = 1;
}

bool ChaseTwiddle::twiddle(int *x, int *y, int *z, std::vector<int> &p)
{
  int i, j, k;

  j = 1;
  while (p[j] <= 0)
    j++;

  if (p[j - 1] == 0) {
    for (int i = j-1; i > 1; --i)
      p[i] = -1;

    p[j] = 0;
    p[1] = 1;

    *x = 0;
    *y = j - 1;
    *z = 0;

  } else {
    if(j > 1)
      p[j - 1] = 0;

    do {
      ++j;
    } while (p[j] > 0);

    i = j;
    k = j - 1;

    while (p[i] == 0)
      p[i++] = -1;

    if (p[i] == -1) {
      p[i] = p[k];

      *x = i - 1;
      *y = k - 1;
      *z = p[k] - 1;

      p[k] = -1;

    } else {
      if (i == p[0]) {
        return true;
      } else {
        p[j] = p[i];

        *x = j - 1;
        *y = i - 1;
        *z = p[i] - 1;

        p[i] = 0;
      }
    }
  }

  return false;
}

} // namespace cgtl
