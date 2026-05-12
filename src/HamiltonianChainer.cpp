#include "HamiltonianChainer.hpp"

#include "ContactMatrix.hpp"
#include "Graph.hpp"

#include <algorithm>
#include <future>
#include <unordered_set>

namespace {

bool is_homozygous(uint32_t id, const BubbleResult& bubbles) {
    return bubbles.phasing_nodes.find(id) == bubbles.phasing_nodes.end();
}

std::vector<uint32_t> collect_subset(const Graph& g,
                                     const BubbleResult& bubbles,
                                     int8_t k) {
    std::vector<uint32_t> v;
    for (uint32_t i = 0; i < g.get_num_nodes(); ++i) {
        const int8_t ph = g.get_phase(i);
        if (ph == k || (ph == -1 && is_homozygous(i, bubbles))) v.push_back(i);
    }
    return v;
}

std::vector<std::vector<uint32_t>>
build_phase_adj(const Graph& g, const std::vector<uint32_t>& subset) {
    std::unordered_set<uint32_t> in_set(subset.begin(), subset.end());
    std::vector<std::vector<uint32_t>> adj(g.get_num_nodes());
    for (uint32_t u : subset) {
        auto range = g.get_neighbors(u);
        for (auto it = range.first; it != range.second; ++it) {
            if (it->source_rev) continue;
            if (in_set.count(it->target_id)) adj[u].push_back(it->target_id);
        }
    }
    return adj;
}

uint32_t edge_weight(uint32_t u, uint32_t v,
                     const std::vector<std::vector<uint32_t>>& adj,
                     const ContactMatrix& contacts) {
    uint32_t w = contacts.get_contact(u, v) * 10u;
    for (uint32_t n : adj[u]) if (n == v) { w += 1; break; }
    return w;
}

// Nearest-neighbor + 2-opt over the phase subgraph.
std::vector<Chain> chain_subset(const Graph& g,
                                const BubbleResult& bubbles,
                                int8_t k,
                                const ContactMatrix& contacts) {
    std::vector<uint32_t> subset = collect_subset(g, bubbles, k);
    if (subset.empty()) return {};

    auto adj = build_phase_adj(g, subset);

    // Seed: in-degree 0 in phase-restricted adj.
    std::unordered_set<uint32_t> in_set(subset.begin(), subset.end());
    std::vector<uint32_t> indeg(g.get_num_nodes(), 0);
    for (uint32_t u : subset)
        for (uint32_t v : adj[u]) ++indeg[v];

    std::vector<char> visited(g.get_num_nodes(), 0);
    std::vector<Chain> chains;

    auto next_seed = [&]() -> uint32_t {
        uint32_t s = UINT32_MAX;
        for (uint32_t u : subset) {
            if (visited[u]) continue;
            if (indeg[u] == 0) return u;
            if (s == UINT32_MAX) s = u;
        }
        return s;
    };

    while (true) {
        uint32_t seed = next_seed();
        if (seed == UINT32_MAX) break;

        std::vector<uint32_t> path;
        uint32_t cur = seed;
        while (cur != UINT32_MAX && !visited[cur]) {
            visited[cur] = 1;
            path.push_back(cur);

            uint32_t best = UINT32_MAX;
            uint32_t best_w = 0;
            for (uint32_t v : adj[cur]) {
                if (visited[v]) continue;
                uint32_t w = edge_weight(cur, v, adj, contacts);
                if (best == UINT32_MAX || w > best_w || (w == best_w && v < best)) {
                    best = v;
                    best_w = w;
                }
            }
            cur = best;
        }

        // 2-opt: swap adjacent segment endpoints if it raises total weight.
        if (path.size() >= 4) {
            auto total = [&](const std::vector<uint32_t>& p) {
                uint64_t s = 0;
                for (size_t i = 0; i + 1 < p.size(); ++i)
                    s += edge_weight(p[i], p[i + 1], adj, contacts);
                return s;
            };
            uint64_t cur_total = total(path);
            bool improved = true;
            while (improved) {
                improved = false;
                for (size_t i = 1; i + 2 < path.size(); ++i) {
                    std::vector<uint32_t> cand = path;
                    std::reverse(cand.begin() + i, cand.begin() + i + 2);
                    uint64_t t = total(cand);
                    if (t > cur_total) {
                        path = std::move(cand);
                        cur_total = t;
                        improved = true;
                        break;
                    }
                }
            }
        }

        Chain c;
        c.phase = k;
        c.nodes = std::move(path);
        c.reverse.assign(c.nodes.size(), false);
        chains.push_back(std::move(c));
    }
    return chains;
}

} // namespace

HamiltonianChainer::HamiltonianChainer(std::shared_ptr<thread_pool::ThreadPool> pool)
    : pool_(std::move(pool)) {}

ChainResult HamiltonianChainer::generate_chain_paths(const Graph& g,
                                                     const BubbleResult& bubbles,
                                                     const ContactMatrix& contacts) const {
    ChainResult res;

    if (pool_) {
        auto f0 = pool_->Submit([&]() { return chain_subset(g, bubbles, 0, contacts); });
        auto f1 = pool_->Submit([&]() { return chain_subset(g, bubbles, 1, contacts); });
        res.phase_0 = f0.get();
        res.phase_1 = f1.get();
    } else {
        res.phase_0 = chain_subset(g, bubbles, 0, contacts);
        res.phase_1 = chain_subset(g, bubbles, 1, contacts);
    }

    for (uint32_t id : bubbles.phasing_nodes) {
        if (g.get_phase(id) == -1) {
            Chain c;
            c.phase = -1;
            c.nodes.push_back(id);
            c.reverse.push_back(false);
            res.unphased.push_back(std::move(c));
        }
    }
    return res;
}
