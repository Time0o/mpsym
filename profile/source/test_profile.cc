#include "partial_perm_inverse_semigroup.h"
#include "partial_perm.h"

using namespace mpsym;

int main()
{
  PartialPermInverseSemigroup M({PartialPerm({1,2},{2,1}),
                                 PartialPerm({1,2},{2,3}),
                                 PartialPerm({3,2},{2,1}),
                                 PartialPerm({1,2,3},{1,2,3})});
                                 // TODO: unity...

//  for (auto const &s : {PartialPerm({1},{2}),
//                        PartialPerm({1},{3}),
//                        PartialPerm({2},{3}),
//                        PartialPerm({1,2},{3,2}),
//                        PartialPerm({2,3},{3,2})})
//    std::cout << s << "/" << (M.contains_element(s) ? "yes" : "no") << std::endl;
//  std::cout << M.contains_element(PartialPerm({2,3},{3,2}));
  M.contains_element(PartialPerm({1,3},{1,2}));

//  std::cout << "===" << std::endl;
//
//  for (auto const &s : {PartialPerm({1,3},{3,1})})
//    std::cout << s << "/" << (M.contains_element(s) ? "yes" : "no") << std::endl;
}
