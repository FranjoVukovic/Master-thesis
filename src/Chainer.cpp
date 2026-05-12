#include "Chainer.hpp"

#include "ContactMatrix.hpp"
#include "Graph.hpp"

#include <algorithm>
#include <unordered_set>

namespace {

bool is_homozygous(uint32_t id, const BubbleResult& bubbles) {
    return bubbles.phasing_nodes.find(id) == bubbles.phasing_nodes.end();
}

// Build the per-haplotype subgraph: nodes with phase==k plus homozygous nodes.
// Edges kept only if both endpoints belong to the subset and target is not the
// alt of any subset node (prevents leaking into opposite haplotype).
std::vector<std::vector<uint32_t>>
build_phase_subgraph(const Graph& g,
                     const BubbleResult& bubbles,
                     int8_t k,
                     std::vector<char>& in_subset) {
    const size_t N = g.get_num_nodes();
    in_subset.assign(N, 0);
    for (uint32_t i = 0; i < N; ++i) {
        const int8_t ph = g.get_phase(i);
        if (ph == k) in_subset[i] = 1;
        else if (ph == -1 && is_homozygous(i, bubbles)) in_subset[i] = 1;
    }

    std::vector<std::vector<uint32_t>> adj(N);
    for (uint32_t u = 0; u < N; ++u) {
        if (!in_subset[u]) continue;
        auto range = g.get_neighbors(u);
        for (auto it = range.first; it != range.second; ++it) {
            if (it->source_rev) continue;            // forward direction only
            uint32_t v = it->target_id;
            if (!in_subset[v]) continue;
            adj[u].push_back(v);
        }
    }
    return adj;
}

uint32_t pick_next(uint32_t u,
                   const std::vector<uint32_t>& cands,
                   const std::vector<char>& visited,
                   const ContactMatrix* contacts) {
    uint32_t best = UINT32_MAX;
    uint32_t best_w = 0;
    for (uint32_t v : cands) {
        if (visited[v]) continue;
        uint32_t w = contacts ? contacts->get_contact(u, v) : 0;
        if (best == UINT32_MAX || w > best_w || (w == best_w && v < best)) {
            best = v;
            best_w = w;
        }
    }
    return best;
}

void chain_haplotype(const Graph& g,
                     const BubbleResult& bubbles,
                     int8_t k,
                     const ContactMatrix* contacts,
                     std::vector<Chain>& out) {
    const size_t N = g.get_num_nodes();
    std::vector<char> in_subset;
    auto adj = build_phase_subgraph(g, bubbles, k, in_subset);

    std::vector<uint32_t> indeg(N, 0);
    for (uint32_t u = 0; u < N; ++u)
        for (uint32_t v : adj[u]) ++indeg[v];

    std::vector<char> visited(N, 0);

    auto walk_from = [&](uint32_t seed) {
        Chain c;
        c.phase = k;
        uint32_t cur = seed;
        while (cur != UINT32_MAX && !visited[cur]) {
            visited[cur] = 1;
            c.nodes.push_back(cur);
            c.reverse.push_back(false);
            cur = pick_next(cur, adj[cur], visited, contacts);
        }
        if (!c.nodes.empty()) out.push_back(std::move(c));
    };

    // Seeds: in-degree 0 first.
    for (uint32_t u = 0; u < N; ++u) {
        if (in_subset[u] && !visited[u] && indeg[u] == 0) walk_from(u);
    }
    // Remaining (cycles or orphans).
    for (uint32_t u = 0; u < N; ++u) {
        if (in_subset[u] && !visited[u]) walk_from(u);
    }
}

} // namespace

ChainResult Chainer::generate_chain_paths(const Graph& g,
                                          const BubbleResult& bubbles,
                                          const ContactMatrix* contacts) const {
    ChainResult res;
    chain_haplotype(g, bubbles, 0, contacts, res.phase_0);
    chain_haplotype(g, bubbles, 1, contacts, res.phase_1);

    // Unphased bubble members: one singleton chain each.
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
