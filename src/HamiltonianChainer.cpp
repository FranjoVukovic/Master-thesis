#include "HamiltonianChainer.hpp"

#include "ContactMatrix.hpp"
#include "Graph.hpp"
#include "HamiltonianPath.hpp"

#include <algorithm>
#include <future>
#include <queue>
#include <unordered_set>

namespace {

bool is_homozygous(uint32_t id, const BubbleResult& bubbles) {
    return bubbles.phasing_nodes.find(id) == bubbles.phasing_nodes.end();
}

// Build per-haplotype subgraph adjacency, forward edges only, restricted to
// {phase==k} ∪ {homozygous}.
std::vector<std::vector<uint32_t>>
build_phase_adj(const Graph& g, const BubbleResult& bubbles, int8_t k,
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
            if (it->source_rev) continue;
            if (in_subset[it->target_id]) adj[u].push_back(it->target_id);
        }
    }
    return adj;
}

// Weakly-connected components over the directed subgraph.
std::vector<std::vector<uint32_t>>
connected_components(const std::vector<std::vector<uint32_t>>& adj,
                     const std::vector<char>& in_subset) {
    const size_t N = adj.size();
    std::vector<std::vector<uint32_t>> undirected(N);
    for (uint32_t u = 0; u < N; ++u) {
        for (uint32_t v : adj[u]) {
            undirected[u].push_back(v);
            undirected[v].push_back(u);
        }
    }
    std::vector<char> seen(N, 0);
    std::vector<std::vector<uint32_t>> comps;
    for (uint32_t s = 0; s < N; ++s) {
        if (!in_subset[s] || seen[s]) continue;
        std::vector<uint32_t> comp;
        std::queue<uint32_t> q;
        q.push(s); seen[s] = 1;
        while (!q.empty()) {
            uint32_t u = q.front(); q.pop();
            comp.push_back(u);
            for (uint32_t v : undirected[u])
                if (!seen[v]) { seen[v] = 1; q.push(v); }
        }
        std::sort(comp.begin(), comp.end());
        comps.push_back(std::move(comp));
    }
    return comps;
}

// Greedy linear walk used as fallback when Hamiltonian is intractable.
std::vector<uint32_t> greedy_walk(const std::vector<uint32_t>& nodes,
                                   const std::vector<std::vector<uint32_t>>& adj,
                                   const ContactMatrix& contacts) {
    std::unordered_set<uint32_t> in_comp(nodes.begin(), nodes.end());
    std::vector<uint32_t> indeg(adj.size(), 0);
    for (uint32_t u : nodes)
        for (uint32_t v : adj[u]) if (in_comp.count(v)) ++indeg[v];

    std::vector<char> visited(adj.size(), 0);
    std::vector<uint32_t> path;

    auto next_seed = [&]() -> uint32_t {
        uint32_t fallback = UINT32_MAX;
        for (uint32_t u : nodes) {
            if (visited[u]) continue;
            if (indeg[u] == 0) return u;
            if (fallback == UINT32_MAX) fallback = u;
        }
        return fallback;
    };

    while (path.size() < nodes.size()) {
        uint32_t seed = next_seed();
        if (seed == UINT32_MAX) break;
        uint32_t cur = seed;
        while (cur != UINT32_MAX && !visited[cur]) {
            visited[cur] = 1;
            path.push_back(cur);
            uint32_t best = UINT32_MAX;
            uint32_t best_w = 0;
            for (uint32_t v : adj[cur]) {
                if (visited[v] || !in_comp.count(v)) continue;
                uint32_t w = contacts.get_contact(cur, v);
                if (best == UINT32_MAX || w > best_w || (w == best_w && v < best)) {
                    best = v; best_w = w;
                }
            }
            cur = best;
        }
    }
    return path;
}

// Solve one connected component via exact Hamiltonian DP if small enough,
// otherwise fall back to greedy. Returns the chain of nodes (in walk order).
std::vector<uint32_t> chain_component(const std::vector<uint32_t>& comp,
                                       const std::vector<std::vector<uint32_t>>& adj,
                                       const ContactMatrix& contacts,
                                       size_t max_iters,
                                       size_t exact_size_cap) {
    if (comp.size() > exact_size_cap) return greedy_walk(comp, adj, contacts);

    std::unordered_set<uint32_t> comp_set(comp.begin(), comp.end());
    std::unordered_set<uint32_t> allowed_starts;
    std::unordered_set<uint32_t> allowed_ends;

    std::vector<uint32_t> indeg(adj.size(), 0);
    std::vector<uint32_t> outdeg(adj.size(), 0);
    for (uint32_t u : comp) {
        for (uint32_t v : adj[u]) if (comp_set.count(v)) { ++outdeg[u]; ++indeg[v]; }
    }
    for (uint32_t u : comp) {
        if (indeg[u]  == 0) allowed_starts.insert(u);
        if (outdeg[u] == 0) allowed_ends.insert(u);
    }

    auto succ = [&](uint32_t u) {
        std::vector<uint32_t> out;
        out.reserve(adj[u].size());
        for (uint32_t v : adj[u]) if (comp_set.count(v)) out.push_back(v);
        return out;
    };

    HamiltonianResult r = find_hamiltonian_path(comp, succ, allowed_starts, allowed_ends, max_iters);
    if (r.solved && !r.path.empty()) return r.path;
    return greedy_walk(comp, adj, contacts);
}

std::vector<Chain> chain_haplotype(const Graph& g,
                                    const BubbleResult& bubbles,
                                    int8_t k,
                                    const ContactMatrix& contacts,
                                    size_t max_iters,
                                    size_t exact_size_cap) {
    std::vector<char> in_subset;
    auto adj = build_phase_adj(g, bubbles, k, in_subset);
    auto comps = connected_components(adj, in_subset);

    std::vector<Chain> chains;
    chains.reserve(comps.size());
    for (const auto& comp : comps) {
        std::vector<uint32_t> walk = chain_component(comp, adj, contacts, max_iters, exact_size_cap);
        if (walk.empty()) continue;
        Chain c;
        c.phase = k;
        c.nodes = std::move(walk);
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
    const size_t max_iters       = 5000;
    const size_t exact_size_cap  = 22;

    ChainResult res;
    if (pool_ && pool_->num_threads() >= 2) {
        auto f0 = pool_->Submit([&]() {
            return chain_haplotype(g, bubbles, 0, contacts, max_iters, exact_size_cap);
        });
        auto f1 = pool_->Submit([&]() {
            return chain_haplotype(g, bubbles, 1, contacts, max_iters, exact_size_cap);
        });
        res.phase_0 = f0.get();
        res.phase_1 = f1.get();
    } else {
        res.phase_0 = chain_haplotype(g, bubbles, 0, contacts, max_iters, exact_size_cap);
        res.phase_1 = chain_haplotype(g, bubbles, 1, contacts, max_iters, exact_size_cap);
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
