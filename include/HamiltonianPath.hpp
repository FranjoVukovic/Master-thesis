#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <unordered_set>
#include <vector>

struct HamiltonianResult {
    bool                  solved = false;   // true if solver finished within max_iters
    std::vector<uint32_t> path;              // a witnessing Hamiltonian path (over target_nodes)
    std::vector<uint32_t> unique_prefix;     // shared prefix across all Hamiltonian solutions
    size_t                iterations = 0;
};

// Find a Hamiltonian path visiting exactly `target_nodes` once each, using the
// successor function `succ_fn` to traverse edges. Iteration-capped bitset DP.
// `allowed_starts` empty => any target node may start; same for `allowed_ends`.
// Hard limit: target_nodes.size() <= 22 (state space ~2^22 * 22). Larger sets
// return solved=false immediately.
HamiltonianResult find_hamiltonian_path(
    const std::vector<uint32_t>& target_nodes,
    const std::function<std::vector<uint32_t>(uint32_t)>& succ_fn,
    const std::unordered_set<uint32_t>& allowed_starts = {},
    const std::unordered_set<uint32_t>& allowed_ends   = {},
    size_t max_iters = 5000);
