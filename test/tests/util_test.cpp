#include <tuple>

#include "gmock/gmock.h"

#include "util.hpp"

#include "test_main.cpp"

using namespace mpsym::internal::util;

class PowTest :
  public testing::TestWithParam<std::tuple<unsigned, unsigned, unsigned>> {};

TEST_P(PowTest, CanObtainIntegerPower)
{
  auto inst = GetParam();
  auto base = std::get<0>(inst);
  auto exp = std::get<1>(inst);
  auto res = std::get<2>(inst);

  EXPECT_EQ(res, pow(base, exp))
    << "Integer power calculated correctly.";
}

INSTANTIATE_TEST_SUITE_P(PowTestInstances, PowTest,
  testing::Values(std::make_tuple(10u, 0u, 1u),
                  std::make_tuple(7u, 1u, 7u),
                  std::make_tuple(2u, 3u, 8u),
                  std::make_tuple(4u, 5u, 1024u),
                  std::make_tuple(3u, 7u, 2187u),
                  std::make_tuple(5u, 3u, 125u)));

class FactorialTest :
  public testing::TestWithParam<std::pair<unsigned, unsigned>> {};

TEST_P(FactorialTest, CanObtainIntegerFactorial)
{
  auto inst = GetParam();
  auto in = std::get<0>(inst);
  auto res = std::get<1>(inst);

  EXPECT_EQ(res, factorial(in))
    << "Integer power calculated correctly.";
}

INSTANTIATE_TEST_SUITE_P(FactorialTestInstances, FactorialTest,
  testing::Values(std::make_pair(0u, 1u),
                  std::make_pair(1u, 1u),
                  std::make_pair(5u, 120u),
                  std::make_pair(7u, 5040u)));
