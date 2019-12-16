#include <vector>

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

namespace cgtl
{

std::vector<PermGroup> PermGroup::wreath_decomposition() const
{
  DBG(DEBUG) << "Finding wreath product decomposition for:";
  DBG(DEBUG) << *this;

  auto blocksystems(BlockSystem::non_trivial(*this));
  auto gens(_bsgs.strong_generators());

  for (BlockSystem const &bs : blocksystems) {
    DBG(TRACE) << "Considering block system: " << bs;
    unsigned d = bs.size();

    PermGroup block_permuter = bs.block_permuter(gens);
    DBG(TRACE) << "Block permuter is:";
    DBG(TRACE) << block_permuter;

    std::vector<PermSet> sigma(d);

    sigma[0] = BlockSystem::block_stabilizers(gens, bs[0]);
    DBG(TRACE) << "Block stabilizer of " << bs[0] << " is: " << sigma[0];

    unsigned tmp = PermGroup(degree(), sigma[0]).order();
    if (_order != util::pow(tmp, d) * block_permuter.order()) {
      DBG(TRACE)
        << "Group order equality not satisfied, skipping block system";
      continue;
    }

    DBG(TRACE) << "Stabilizers of remaining blocks are:";
    for (unsigned i = 1u; i < d; ++i) {
      sigma[i] = BlockSystem::block_stabilizers(gens, bs[i]);
      DBG(TRACE) << bs[i] << " => " << sigma[i];
    }

    // try to find monomorphism from block permuter to group using heuristic
    PermSet block_permuter_generator_image;

    std::vector<unsigned> tmp_perm(degree());
    for (Perm const &gen : block_permuter.bsgs().strong_generators()) {
      for (unsigned i = 0u; i < d; ++i) {
        auto block(bs[i]);

        for (auto j = 0u; j < block.size(); ++j)
          tmp_perm[block[j] - 1u] = bs[gen[i + 1u] - 1u][j];
      }

      block_permuter_generator_image.emplace(tmp_perm);
    }

    DBG(TRACE) << "Heuristic monomorphism image is:";
    DBG(TRACE) << block_permuter_generator_image;

    bool found_monomorphism = true;

    auto classes(bs.classes());

    PermSet block_permuter_generator_reconstruction;

    tmp_perm.resize(d);
    for (Perm const &gen : block_permuter_generator_image) {
      for (unsigned i = 0u; i < d; ++i) {
        unsigned x = bs[i][0];
        unsigned y = gen[x];

        tmp_perm[i] = classes[y - 1u];
      }

      Perm reconstructed_gen(tmp_perm);

      if (!block_permuter.contains_element(reconstructed_gen)) {
        found_monomorphism = false;
        break;
      }

      block_permuter_generator_reconstruction.insert(reconstructed_gen);
    }

    DBG(TRACE) << "Block permuter reconstruction yields generators:";
    DBG(TRACE) << block_permuter_generator_reconstruction;

    if (found_monomorphism) {
      PermGroup block_permuter_reconstructed(
        d, block_permuter_generator_reconstruction);

      if (block_permuter_reconstructed.order() != block_permuter.order())
        found_monomorphism = false;
    }

    if (!found_monomorphism) {
      DBG(WARN)
        << "Wreath decomposition exists but was not found by heuristic";

      break;
    }

    std::vector<PermGroup> res(d + 1u);

    res[0] = PermGroup(degree(), block_permuter_generator_image);
    for (unsigned i = 0u; i < d; ++i)
      res[i + 1u] = PermGroup(degree(), sigma[i]);

    DBG(TRACE)
      << "==> Found wreath product decomposition, listing generators:";
#ifndef NDEBUG
    for (PermGroup const &pg : res)
      DBG(TRACE) << pg.bsgs().strong_generators();
#endif

    return res;
  }

  DBG(TRACE) << "==> No wreath product decomposition found";
  return std::vector<PermGroup>();
}

} // namespace cgtl
