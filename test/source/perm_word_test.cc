#include <unordered_set>
#include <vector>

#include "gmock/gmock.h"

#include "perm.h"
#include "test_utility.h"

#include "test_main.cc"

using cgtl::Perm;
using cgtl::PermWord;

using testing::UnorderedElementsAreArray;

TEST(PermWordTest, CanConstructPermWord)
{
  PermWord perm;
  EXPECT_TRUE(perm_word_equal({1}, perm))
    << "Default construction produces identity permutation word.";

  PermWord perm_id(5);
  EXPECT_TRUE(perm_word_equal({1, 2, 3, 4, 5}, perm_id))
    << "Identity construction produces identity permutation word.";

  PermWord perm_explicit({1, 3, 4, 5, 2});
  EXPECT_TRUE(perm_word_equal({1, 3, 4, 5, 2}, perm_explicit))
    << "Explicit construction produces correct permutation word.";

  PermWord perm_empty_cycle(6, {});
  EXPECT_TRUE(perm_word_equal({1, 2, 3, 4, 5, 6}, perm_empty_cycle))
    << "No-cycles construction produces correct permutation word.";

  PermWord perm_single_cycle(6, {{3, 2, 5}});
  EXPECT_TRUE(perm_word_equal({1, 5, 2, 4, 3, 6}, perm_single_cycle))
    << "Single-cycle construction produces correct permutation word.";

  PermWord perm_multi_cycles(6, {{6, 2, 4}, {2, 5, 4}, {3, 2, 5}});
  EXPECT_TRUE(perm_word_equal({1, 5, 2, 6, 4, 3}, perm_multi_cycles))
    << "Multi-cycle construction produces correct permutation word.";

  Perm simple_perm_multi_cycles(6, {{6, 2, 4}, {2, 5, 4}, {3, 2, 5}});
  PermWord perm_word_multi_cycles(simple_perm_multi_cycles);
  EXPECT_TRUE(perm_word_equal({1, 5, 2, 6, 4, 3}, perm_word_multi_cycles))
    << "Can convert simple permutation to permutation word.";

  EXPECT_TRUE(perm_equal({1, 5, 2, 6, 4, 3}, perm_word_multi_cycles.perm()))
    << "Can obtain simple permutation from permutation word.";
}

TEST(PermWordTest, CanInvertPermWord)
{
  PermWord perm({3, 2, 4, 1});

  EXPECT_TRUE(perm_word_equal({4, 2, 1, 3}, ~perm))
    << "Inverting permutation word works.";

  PermWord perm1({3, 2, 4, 1});
  PermWord perm2({1, 4, 3, 2});
  EXPECT_TRUE(perm_word_equal({4, 3, 1, 2}, ~(perm1 * perm2)))
    << "Inverting permutation words of even length works.";

  PermWord perm3({3, 2, 4, 1});
  PermWord perm4({2, 1, 3, 4});
  PermWord perm5({4, 3, 2, 1});

  EXPECT_TRUE(perm_word_equal({3, 1, 4, 2}, ~(perm3 * perm4 * perm5)))
    << "Inverting permutation words of uneven length works.";
}

TEST(PermWordTest, CanMultiplyPermWords)
{
  PermWord perm0(7, {{1, 2, 4}});
  perm0 *= PermWord(7, {{4, 5}});

  EXPECT_TRUE(perm_word_equal({2, 5, 3, 1, 4, 6, 7}, perm0))
    << "Multiplying plus assigning permutation words produces correct result.";

  PermWord perm1(6, {{2, 5, 4}});
  PermWord perm2(6, {{3, 2, 5}});

  PermWord perm_mult1 = perm1 * perm2;
  EXPECT_TRUE(perm_word_equal({1, 3, 2, 5, 4, 6}, perm_mult1))
    << "Multiplying permutation words produces correct result.";
}

TEST(PermWordTest, PermWordStringRepresentation)
{
  PermWord perm1({2, 3, 1, 5, 4});
  std::stringstream ss1;
  ss1 << perm1;

  EXPECT_EQ("(1 2 3)(4 5)", ss1.str())
    << "Correct permutation word string representation.";

  PermWord perm2({1, 5, 3, 6, 2, 7, 4, 8});
  std::stringstream ss2;
  ss2 << perm2;

  EXPECT_EQ("(2 5)(4 6 7)", ss2.str())
    << "Permutation word string representation ignores single element cycles.";

  PermWord perm3({1, 2, 3});
  std::stringstream ss3;
  ss3 << perm3;

  EXPECT_EQ("()", ss3.str())
    << "Identity permutation string representation correct.";
}

TEST(PermWordTest, CanHashPermWord)
{
  std::vector<PermWord> perm_words = {
    PermWord(5, {{1, 2, 3}}),
    PermWord(5, {{2, 3}, {4, 5}}),
    PermWord(5, {{1, 2, 3, 4}}),
    PermWord(5, {{1, 2}}),
    PermWord(5, {{1, 2, 3}, {4, 5}})
  };

  std::unordered_set<PermWord> perm_word_set;

  unsigned const repetitions = 10u;
  for (unsigned i = 0u; i < repetitions; ++i) {
    for (PermWord const &pw : perm_words)
      perm_word_set.insert(pw);
  }

  ASSERT_EQ(perm_words.size(), perm_word_set.size())
    << "Hashed permutation word set has correct size.";

  std::vector<PermWord>hashed_perm_words(perm_word_set.begin(),
                                         perm_word_set.end());

  EXPECT_THAT(hashed_perm_words, UnorderedElementsAreArray(perm_words))
    << "Hashed permutation word set has correct elements.";
}
