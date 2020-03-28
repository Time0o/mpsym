#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#include "block_system.h"
#include "dbg.h"
#include "perm.h"
#include "perm_group.h"
#include "perm_set.h"
#include "util.h"

/**
 * @file perm_group_wreath_decomp.h
 * @brief Implements `PermGroup::wreath_decomposition`.
 *
 * @author Timo Nicolai
 */

namespace mpsym
{

std::vector<PermGroup> PermGroup::wreath_decomposition() const
{
  DBG(DEBUG) << "Finding wreath product decomposition for:";
  DBG(DEBUG) << *this;

  for (BlockSystem const &block_system : BlockSystem::non_trivial(*this)) {
    DBG(TRACE) << "Considering block system:";
    DBG(TRACE) << block_system;

    // determine block permuter subgroup
    PermGroup block_permuter(block_system.size(),
                             block_system.block_permuter(generators()));

    DBG(TRACE) << "Block permuter is:";
    DBG(TRACE) << block_permuter;

    // determine block stabilizer subgroups
    auto stabilizers(wreath_decomp_find_stabilizers(
      block_system, block_permuter));

    if (stabilizers.empty())
      continue;

    // check if a monomorphism can be found heuristically
    auto block_permuter_image(wreath_decomp_construct_block_permuter_image(
      block_system, block_permuter));

    bool found_monomorphism(wreath_decomp_reconstruct_block_permuter(
      block_system, block_permuter, block_permuter_image));

    if (!found_monomorphism)
      break;

    // construct the wreath decomposition
    std::vector<PermGroup> decomposition(block_system.size() + 1u);

    decomposition[0] = PermGroup(degree(), block_permuter_image);
    for (unsigned i = 0u; i < block_system.size(); ++i)
      decomposition[i + 1u] = stabilizers[i];

    DBG(DEBUG) << "==> Found wreath product decomposition:";
#ifndef NDEBUG
    for (PermGroup const &pg : decomposition)
      DBG(DEBUG) << pg;
#endif

    return decomposition;
  }

  DBG(DEBUG) << "==> No wreath product decomposition found";
  return {};
}

std::vector<PermGroup> PermGroup::wreath_decomp_find_stabilizers(
  BlockSystem const &block_system,
  PermGroup const &block_permuter) const
{
  using boost::multiprecision::pow;

  std::vector<PermGroup> stabilizers(block_system.size());

  auto create_stabilizer = [&](unsigned i) {
    auto block(block_system[i]);

    // find stabilizer subgroup generators
    auto stabilizer_generators(
      BlockSystem::block_stabilizers(generators(), block));

    // restrict stabilizer subgroup generators to block
    PermSet stabilizer_generators_restricted;
    for (Perm const &gen : stabilizer_generators) {
      stabilizer_generators_restricted.insert(
        gen.restricted(block.begin(), block.end()));
    }

    // construct stabilizer subgroup
    stabilizers[i] = PermGroup(degree(), stabilizer_generators_restricted);

    DBG(TRACE) << "Block stabilizer of " << block_system[i] << ":";
    DBG(TRACE) << stabilizers[i];
  };

  // determine stabilizer subgroup of first block
  create_stabilizer(0);

  // skip blocksystem if order equality not fullfilled
  BSGS::order_type expected_order = pow(
    stabilizers[0].order(), block_system.size()) * block_permuter.order();

  if (_order != expected_order) {
    DBG(TRACE) << "Group order equality not satisfied";
    return {};
  }

  // determine stabilizers subgroups of remaining blocks
  for (unsigned i = 1u; i < block_system.size(); ++i)
    create_stabilizer(i);

  return stabilizers;
}

PermSet PermGroup::wreath_decomp_construct_block_permuter_image(
  BlockSystem const &block_system,
  PermGroup const &block_permuter) const
{
  PermSet block_permuter_image;

  for (Perm const &gen : block_permuter.generators()) {
    std::vector<unsigned> perm(degree());

    for (unsigned i = 0u; i < block_system.size(); ++i) {
      auto block(block_system[i]);

      for (auto j = 0u; j < block.size(); ++j)
        perm[block[j] - 1u] = block_system[gen[i + 1u] - 1u][j];
    }

    block_permuter_image.insert(Perm(perm));
  }

  DBG(TRACE) << "Heuristic monomorphism image generators:";
  DBG(TRACE) << block_permuter_image;

  return block_permuter_image;
}

bool PermGroup::wreath_decomp_reconstruct_block_permuter(
  BlockSystem const &block_system,
  PermGroup const &block_permuter,
  PermSet const &block_permuter_image) const
{
  bool found_monomorphism = true;

  PermSet block_permuter_reconstruction;

  for (Perm const &gen : block_permuter_image) {
    std::vector<unsigned> perm(block_system.size());

    for (unsigned i = 0u; i < block_system.size(); ++i)
      perm[i] = block_system.block_index(gen[block_system[i][0]]) + 1u;

    Perm reconstructed_gen(perm);

    if (!block_permuter.contains_element(reconstructed_gen)) {
      found_monomorphism = false;
      break;
    }

    block_permuter_reconstruction.insert(reconstructed_gen);
  }

  DBG(TRACE) << "Block permuter reconstruction yields generators:";
  DBG(TRACE) << block_permuter_reconstruction;

  if (found_monomorphism) {
    if (PermGroup(block_system.size(), block_permuter_reconstruction).order()
        != block_permuter.order())
      found_monomorphism = false;
  }

  if (!found_monomorphism) {
    DBG(WARN) << "Wreath decomposition exists but was not found by heuristic";
  }

  return found_monomorphism;
}

} // namespace mpsym
