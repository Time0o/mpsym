#include <iostream>

#include "bsgs.h"
#include "perm_group.h"
#include "perm_set.h"

using namespace cgtl;

int main()
{
  PermGroup pg1(12, {Perm(12, {{1,2}}), Perm(12, {{1,2,3,4,5,6,7,8,9,10,11,12}})});

  std::cout << pg1.order() << std::endl;

  PermGroup pg2(12,
    {
      Perm(12, {{1,2,3}}),
      Perm(12, {{1,2,4}}),
      Perm(12, {{1,2,5}}),
      Perm(12, {{1,2,6}}),
      Perm(12, {{1,2,7}}),
      Perm(12, {{1,2,8}}),
      Perm(12, {{1,2,9}}),
      Perm(12, {{1,2,10}}),
      Perm(12, {{1,2,11}}),
      Perm(12, {{1,2,12}})
    }
  );

  std::cout << pg2.order() << std::endl;
}
