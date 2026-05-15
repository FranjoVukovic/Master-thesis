#include "HamiltonianPath.hpp"

#include <algorithm>
#include <unordered_map>

namespace {

// Reconstruct a Hamiltonian path by retracing predecessor pointers in the DP.
// `pred[(mask, end)]` = (prev_mask, prev_end) or sentinel for a start state.
struct State {
    uint32_t mask;
    uint32_t end_idx;
    bool operator==(const State& o) const { return mask == o.mask && end_idx == o.end_idx; }
};
struct StateHash {
    size_t operator()(const State& s) const {
        return (static_cast<size_t>(s.mask) << 8) ^ s.end_idx;
    }
};

std::vector<uint32_t> reconstruct(
    State final_state,
    const std::unordered_map<State, State, StateHash>& pred,
    const std::vector<uint32_t>& node_of_idx) {

    std::vector<uint32_t> rev;
    State cur = final_state;
    while (true) {
        rev.push_back(node_of_idx[cur.end_idx]);
        auto it = pred.find(cur);
        if (it == pred.end()) break;  // start state (no predecessor recorded)
        cur = it->second;
    }
    std::reverse(rev.begin(), rev.end());
    return rev;
}

} // namespace

HamiltonianResult find_hamiltonian_path(
    const std::vector<uint32_t>& target_nodes,
    const std::function<std::vector<uint32_t>(uint32_t)>& succ_fn,
    const std::unordered_set<uint32_t>& allowed_starts,
    const std::unordered_set<uint32_t>& allowed_ends,
    size_t max_iters) {

    HamiltonianResult res;
    const size_t n = target_nodes.size();
    if (n == 0) { res.solved = true; return res; }
    if (n > 22) return res;  // too big for exact DP
    if (n == 1) {
        const uint32_t v = target_nodes[0];
        const bool ok_start = allowed_starts.empty() || allowed_starts.count(v);
        const bool ok_end   = allowed_ends.empty()   || allowed_ends.count(v);
        if (ok_start && ok_end) {
            res.solved = true;
            res.path = {v};
            res.unique_prefix = {v};
        }
        return res;
    }

    // Build local index for bitmask.
    std::unordered_map<uint32_t, uint32_t> idx_of;
    idx_of.reserve(n * 2);
    for (uint32_t i = 0; i < n; ++i) idx_of[target_nodes[i]] = i;

    // Adjacency restricted to target_nodes, indexed.
    std::vector<std::vector<uint32_t>> adj(n);
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t v : succ_fn(target_nodes[i])) {
            auto it = idx_of.find(v);
            if (it != idx_of.end()) adj[i].push_back(it->second);
        }
    }

    const uint32_t full_mask = (n == 32) ? 0xFFFFFFFFu
                                         : ((static_cast<uint32_t>(1) << n) - 1u);

    // DP: forward expansion via BFS over (mask, end) states.
    // `visited` tracks reachability; `pred` only records non-start states.
    std::unordered_set<State, StateHash> visited;
    std::unordered_map<State, State, StateHash> pred;
    std::vector<State> frontier;
    frontier.reserve(64);
    for (uint32_t i = 0; i < n; ++i) {
        if (!allowed_starts.empty() && !allowed_starts.count(target_nodes[i])) continue;
        State s{static_cast<uint32_t>(1u << i), i};
        if (visited.insert(s).second) frontier.push_back(s);
    }

    size_t iters = 0;
    std::vector<State> next_frontier;
    while (!frontier.empty() && iters < max_iters) {
        next_frontier.clear();
        for (const State& s : frontier) {
            for (uint32_t nbr : adj[s.end_idx]) {
                const uint32_t bit = static_cast<uint32_t>(1u) << nbr;
                if (s.mask & bit) continue;  // already visited in this path
                State t{s.mask | bit, nbr};
                if (visited.insert(t).second) {
                    pred[t] = s;
                    next_frontier.push_back(t);
                    ++iters;
                    if (iters >= max_iters) break;
                }
            }
            if (iters >= max_iters) break;
        }
        frontier.swap(next_frontier);
    }
    res.iterations = iters;

    // Collect all full-mask states whose end is in allowed_ends.
    std::vector<State> winners;
    for (const State& s : visited) {
        if (s.mask != full_mask) continue;
        if (!allowed_ends.empty() && !allowed_ends.count(target_nodes[s.end_idx])) continue;
        winners.push_back(s);
    }
    if (winners.empty()) {
        // No Hamiltonian path — either intractable or none exists.
        // Treat as solved=true only if iteration budget wasn't exhausted; otherwise unsolved.
        res.solved = (iters < max_iters);
        return res;
    }

    res.solved = true;
    res.path = reconstruct(winners.front(), pred, target_nodes);

    // Unique prefix: longest common prefix across all winning paths.
    std::vector<uint32_t> common = res.path;
    for (size_t i = 1; i < winners.size(); ++i) {
        auto p = reconstruct(winners[i], pred, target_nodes);
        size_t lcp = 0;
        const size_t lim = std::min(common.size(), p.size());
        while (lcp < lim && common[lcp] == p[lcp]) ++lcp;
        common.resize(lcp);
        if (common.empty()) break;
    }
    res.unique_prefix = std::move(common);
    return res;
}
