#ifndef _GUARD_CHASE_TWIDDLE
#define _GUARD_CHASE_TWIDDLE

#include <vector>

namespace cgtl
{

class ChaseTwiddle
{
public:
  ChaseTwiddle(unsigned n, unsigned k)
  : _n(n),
    _k(k)
  {}

  template<typename T, typename FUNC>
  void foreach(T const &a_, FUNC &&f)
  {
    T a(a_);
    f(a);

    std::vector<int> p(_n + 2);

    init(_n, _k, p);

    int x, y, z;
    while (!twiddle(&x, &y, &z, p)) {
      a[z] = x;
      f(a);
    }
  }

private:
  void init(int n, int k, std::vector<int> &p);
  bool twiddle(int *x, int *y, int *z, std::vector<int> &p);

  unsigned _n, _k;
  std::vector<int> _p;
};

} // namespace cgtl

#endif // _GUARD_CHASE_TWIDDLE
