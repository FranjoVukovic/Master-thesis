#include "Phaser.hpp"

#include "ContactMatrix.hpp"
#include "Graph.hpp"

#include <algorithm>
#include <cstdint>
#include <future>
#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// ============================================================
// Trio-constrained initialization
// ============================================================

size_t apply_trio_constraints(Graph& graph,
                               const std::unordered_set<uint32_t>& phasing_nodes,
                               const std::vector<TrioScores>&      scores,
                               const PhaserConfig&                 cfg) {
    size_t locked = 0;
    for (uint32_t id : phasing_nodes) {
        if (id >= scores.size()) continue;
        const TrioScores& s = scores[id];
        const uint64_t total = static_cast<uint64_t>(s.pat_count) + s.mat_count;
        if (total < cfg.min_total_kmers) continue;

        const double pat_frac = static_cast<double>(s.pat_count) / static_cast<double>(total);
        const double mat_frac = static_cast<double>(s.mat_count) / static_cast<double>(total);

        if (pat_frac >= cfg.lock_ratio) {
            graph.lock_phase(id, 0);
            ++locked;
        } else if (mat_frac >= cfg.lock_ratio) {
            graph.lock_phase(id, 1);
            ++locked;
        }
    }
    return locked;
}

// ============================================================
// Monte Carlo phasing
// ============================================================

namespace {

constexpr size_t NONE_META   = static_cast<size_t>(-1);
constexpr size_t LOCKED_META = static_cast<size_t>(-2);

struct MetaBubble {
    uint32_t              a = 0;
    uint32_t              b = 0;
    std::vector<uint32_t> member_nodes_a;
    std::vector<uint32_t> member_nodes_b;
    std::vector<size_t>   orig_indices;
};

class UnionFind {
public:
    explicit UnionFind(size_t n) : parent_(n), rank_(n, 0) {
        std::iota(parent_.begin(), parent_.end(), size_t{0});
    }
    size_t find(size_t x) {
        while (parent_[x] != x) {
            parent_[x] = parent_[parent_[x]];
            x = parent_[x];
        }
        return x;
    }
    void unite(size_t x, size_t y) {
        size_t rx = find(x), ry = find(y);
        if (rx == ry) return;
        if (rank_[rx] < rank_[ry]) std::swap(rx, ry);
        parent_[ry] = rx;
        if (rank_[rx] == rank_[ry]) ++rank_[rx];
    }
private:
    std::vector<size_t> parent_;
    std::vector<size_t> rank_;
};

void build_node_maps(const std::vector<MetaBubble>&      meta_list,
                     const Graph&                        graph,
                     const std::unordered_set<uint32_t>& phasing_nodes,
                     std::vector<size_t>&                node_to_bubble_idx,
                     std::vector<int8_t>&                node_side) {
    std::fill(node_to_bubble_idx.begin(), node_to_bubble_idx.end(), NONE_META);
    std::fill(node_side.begin(), node_side.end(), int8_t{0});
    for (size_t i = 0; i < meta_list.size(); ++i) {
        for (uint32_t u : meta_list[i].member_nodes_a) {
            node_to_bubble_idx[u] = i;
            node_side[u]          = +1;
        }
        for (uint32_t u : meta_list[i].member_nodes_b) {
            node_to_bubble_idx[u] = i;
            node_side[u]          = -1;
        }
    }
    // Trio-locked phasing nodes contribute fixed signs to ΔS.
    for (uint32_t id : phasing_nodes) {
        if (id >= node_to_bubble_idx.size()) continue;
        if (node_to_bubble_idx[id] != NONE_META) continue;
        if (graph.is_phase_locked(id)) {
            node_to_bubble_idx[id] = LOCKED_META;
            node_side[id] = (graph.get_phase(id) == 0) ? int8_t{+1} : int8_t{-1};
        }
    }
}

std::vector<std::vector<size_t>> build_bubble_adj(
        const std::vector<MetaBubble>& meta_list,
        const ContactMatrix&           contacts,
        const std::vector<size_t>&     node_to_bubble_idx) {
    std::vector<std::vector<size_t>> adj(meta_list.size());
    for (size_t i = 0; i < meta_list.size(); ++i) {
        std::unordered_set<size_t> seen;
        auto walk = [&](uint32_t u) {
            const auto end = contacts.row_end(u);
            for (auto pos = contacts.row_begin(u); pos != end; ++pos) {
                uint32_t n = contacts.col_at(pos);
                size_t   j = node_to_bubble_idx[n];
                if (j >= meta_list.size()) continue;   // NONE or LOCKED
                if (j == i) continue;
                seen.insert(j);
            }
        };
        for (uint32_t u : meta_list[i].member_nodes_a) walk(u);
        for (uint32_t u : meta_list[i].member_nodes_b) walk(u);
        adj[i].assign(seen.begin(), seen.end());
        std::sort(adj[i].begin(), adj[i].end());
    }
    return adj;
}

double delta_for_flip(size_t                          meta_idx,
                      int8_t                          sa,
                      const std::vector<MetaBubble>&  meta_list,
                      const std::vector<int8_t>&      signs,
                      const ContactMatrix&            contacts,
                      const std::vector<size_t>&      node_to_bubble_idx,
                      const std::vector<int8_t>&      node_side) {
    int64_t sum = 0;
    const MetaBubble& meta = meta_list[meta_idx];
    auto contribute = [&](uint32_t u, int8_t su) {
        const auto end = contacts.row_end(u);
        for (auto pos = contacts.row_begin(u); pos != end; ++pos) {
            const uint32_t n = contacts.col_at(pos);
            const uint32_t w = contacts.value_at(pos);
            const size_t   j = node_to_bubble_idx[n];
            if (j == NONE_META) continue;
            if (j == meta_idx)  continue;
            int8_t sj;
            if (j == LOCKED_META) sj = node_side[n];
            else                  sj = static_cast<int8_t>(signs[j] * node_side[n]);
            sum += static_cast<int64_t>(w) * su * sj;
        }
    };
    for (uint32_t u : meta.member_nodes_a) contribute(u, sa);
    for (uint32_t u : meta.member_nodes_b) contribute(u, static_cast<int8_t>(-sa));
    return -2.0 * static_cast<double>(sum);
}

std::vector<int8_t> run_sample(uint64_t                       seed,
                                int                            iterations,
                                const std::vector<MetaBubble>& meta_list,
                                const ContactMatrix&           contacts,
                                const std::vector<size_t>&     node_to_bubble_idx,
                                const std::vector<int8_t>&     node_side) {
    std::mt19937_64 rng(seed);
    std::vector<int8_t> signs(meta_list.size());
    for (size_t i = 0; i < signs.size(); ++i) {
        signs[i] = (rng() & 1ULL) ? int8_t{+1} : int8_t{-1};
    }
    for (int it = 0; it < iterations; ++it) {
        for (size_t i = 0; i < meta_list.size(); ++i) {
            const double d = delta_for_flip(i, signs[i], meta_list, signs,
                                            contacts, node_to_bubble_idx, node_side);
            if (d > 0) signs[i] = static_cast<int8_t>(-signs[i]);
        }
    }
    return signs;
}

} // namespace

void monte_carlo_phase(Graph&                                   graph,
                       const ContactMatrix&                     contacts,
                       const BubbleResult&                      bubbles,
                       const PhaserConfig&                      cfg,
                       std::shared_ptr<thread_pool::ThreadPool> pool) {
    // 1. FREE_PAIRS: canonical (a<b), neither node trio-locked.
    std::vector<std::pair<uint32_t, uint32_t>> original_pairs;
    original_pairs.reserve(bubbles.alt_map.size() / 2);
    for (const auto& kv : bubbles.alt_map) {
        const uint32_t a = kv.first, b = kv.second;
        if (a >= b) continue;
        if (graph.is_phase_locked(a) || graph.is_phase_locked(b)) continue;
        original_pairs.emplace_back(a, b);
    }
    if (original_pairs.empty()) return;
    std::sort(original_pairs.begin(), original_pairs.end());

    const size_t N_orig  = original_pairs.size();
    const size_t N_nodes = graph.get_num_nodes();

    // 2. Initial meta_list (one meta per original bubble).
    std::vector<MetaBubble> meta_list(N_orig);
    std::vector<int8_t>     orig_sign_in_meta(N_orig, +1);
    std::vector<size_t>     meta_of_orig(N_orig);
    for (size_t i = 0; i < N_orig; ++i) {
        meta_list[i].a              = original_pairs[i].first;
        meta_list[i].b              = original_pairs[i].second;
        meta_list[i].member_nodes_a = {original_pairs[i].first};
        meta_list[i].member_nodes_b = {original_pairs[i].second};
        meta_list[i].orig_indices   = {i};
        meta_of_orig[i] = i;
    }

    const size_t effective_samples = std::max<size_t>(1, pool->num_threads());

    std::vector<size_t> node_to_bubble_idx(N_nodes, NONE_META);
    std::vector<int8_t> node_side(N_nodes, 0);

    // 3. Rounds: parallel sample → agreement merge → rebuild metas.
    for (int round_idx = 0; round_idx < cfg.n_rounds; ++round_idx) {
        if (meta_list.size() <= 1) break;

        build_node_maps(meta_list, graph, bubbles.phasing_nodes,
                        node_to_bubble_idx, node_side);
        auto bubble_adj = build_bubble_adj(meta_list, contacts, node_to_bubble_idx);

        std::vector<std::future<std::vector<int8_t>>> futs;
        futs.reserve(effective_samples);
        for (size_t s = 0; s < effective_samples; ++s) {
            const uint64_t seed = static_cast<uint64_t>(cfg.rng_seed) ^
                                  (static_cast<uint64_t>(round_idx) * 1000ULL +
                                   static_cast<uint64_t>(s));
            const int iters = cfg.core_iterations;
            futs.emplace_back(pool->Submit(
                [seed, iters,
                 &meta_list, &contacts, &node_to_bubble_idx, &node_side]() {
                    return run_sample(seed, iters, meta_list, contacts,
                                      node_to_bubble_idx, node_side);
                }));
        }
        std::vector<std::vector<int8_t>> all_signs(effective_samples);
        for (size_t s = 0; s < effective_samples; ++s) {
            all_signs[s] = futs[s].get();
        }

        // Agreement merge over contact-adjacent meta pairs only.
        UnionFind uf(meta_list.size());
        for (size_t i = 0; i < meta_list.size(); ++i) {
            for (size_t j : bubble_adj[i]) {
                if (j <= i) continue;
                size_t n_par = 0, n_anti = 0;
                for (size_t s = 0; s < effective_samples; ++s) {
                    if (all_signs[s][i] == all_signs[s][j]) ++n_par;
                    else                                    ++n_anti;
                }
                const double agreement =
                    static_cast<double>(std::max(n_par, n_anti)) /
                    static_cast<double>(effective_samples);
                if (agreement >= cfg.merge_threshold) uf.unite(i, j);
            }
        }

        // Group metas by component root, sorted for deterministic order.
        std::unordered_map<size_t, std::vector<size_t>> comp_members;
        for (size_t i = 0; i < meta_list.size(); ++i) {
            comp_members[uf.find(i)].push_back(i);
        }
        std::vector<std::vector<size_t>> sorted_comps;
        sorted_comps.reserve(comp_members.size());
        for (auto& kv : comp_members) {
            std::sort(kv.second.begin(), kv.second.end());
            sorted_comps.push_back(std::move(kv.second));
        }
        std::sort(sorted_comps.begin(), sorted_comps.end(),
                  [](const std::vector<size_t>& x, const std::vector<size_t>& y) {
                      return x.front() < y.front();
                  });

        // Rebuild meta list. For each component, rep = smallest meta index;
        // non-rep members get a relative sign by majority vote, swapping their
        // a/b sides if rel == -1.
        std::vector<MetaBubble> new_meta_list;
        new_meta_list.reserve(sorted_comps.size());
        std::vector<int8_t> new_orig_sign(N_orig, +1);
        std::vector<size_t> new_meta_of_orig(N_orig);

        for (const std::vector<size_t>& members : sorted_comps) {
            const size_t rep = members.front();
            const size_t new_idx = new_meta_list.size();

            MetaBubble m;
            m.member_nodes_a = meta_list[rep].member_nodes_a;
            m.member_nodes_b = meta_list[rep].member_nodes_b;
            m.orig_indices   = meta_list[rep].orig_indices;
            for (size_t o : meta_list[rep].orig_indices) {
                new_orig_sign[o]    = orig_sign_in_meta[o];
                new_meta_of_orig[o] = new_idx;
            }

            for (size_t k = 1; k < members.size(); ++k) {
                const size_t mb = members[k];
                int votes_pos = 0, votes_neg = 0;
                for (size_t s = 0; s < effective_samples; ++s) {
                    if (all_signs[s][mb] == all_signs[s][rep]) ++votes_pos;
                    else                                       ++votes_neg;
                }
                const int8_t rel = (votes_pos >= votes_neg) ? int8_t{+1} : int8_t{-1};

                if (rel == +1) {
                    m.member_nodes_a.insert(m.member_nodes_a.end(),
                        meta_list[mb].member_nodes_a.begin(),
                        meta_list[mb].member_nodes_a.end());
                    m.member_nodes_b.insert(m.member_nodes_b.end(),
                        meta_list[mb].member_nodes_b.begin(),
                        meta_list[mb].member_nodes_b.end());
                } else {
                    m.member_nodes_a.insert(m.member_nodes_a.end(),
                        meta_list[mb].member_nodes_b.begin(),
                        meta_list[mb].member_nodes_b.end());
                    m.member_nodes_b.insert(m.member_nodes_b.end(),
                        meta_list[mb].member_nodes_a.begin(),
                        meta_list[mb].member_nodes_a.end());
                }
                m.orig_indices.insert(m.orig_indices.end(),
                    meta_list[mb].orig_indices.begin(),
                    meta_list[mb].orig_indices.end());
                for (size_t o : meta_list[mb].orig_indices) {
                    new_orig_sign[o]    = static_cast<int8_t>(orig_sign_in_meta[o] * rel);
                    new_meta_of_orig[o] = new_idx;
                }
            }
            m.a = m.member_nodes_a.front();
            m.b = m.member_nodes_b.front();
            new_meta_list.push_back(std::move(m));
        }

        meta_list         = std::move(new_meta_list);
        orig_sign_in_meta = std::move(new_orig_sign);
        meta_of_orig      = std::move(new_meta_of_orig);
    }

    // 4. Final greedy pass on reduced graph (single-threaded, 3x iterations).
    build_node_maps(meta_list, graph, bubbles.phasing_nodes,
                    node_to_bubble_idx, node_side);
    const uint64_t final_seed = static_cast<uint64_t>(cfg.rng_seed) ^
                                (static_cast<uint64_t>(cfg.n_rounds) * 1000ULL);
    const std::vector<int8_t> final_signs =
        run_sample(final_seed, 3 * cfg.core_iterations, meta_list, contacts,
                   node_to_bubble_idx, node_side);

    // 5. Writeback: combine per-orig sign with its meta's final sign.
    for (size_t o = 0; o < N_orig; ++o) {
        const int8_t combined = static_cast<int8_t>(
            orig_sign_in_meta[o] * final_signs[meta_of_orig[o]]);
        const uint32_t a = original_pairs[o].first;
        const uint32_t b = original_pairs[o].second;
        const int8_t pa = (combined == +1) ? int8_t{0} : int8_t{1};
        const int8_t pb = (combined == +1) ? int8_t{1} : int8_t{0};
        if (!graph.is_phase_locked(a)) graph.lock_phase(a, pa);
        if (!graph.is_phase_locked(b)) graph.lock_phase(b, pb);
    }

    // 6. Propagate to interior multi-branch members: each alt_map key's phase
    //    is the opposite of its alt target. Iterate to convergence (≤2 sweeps
    //    in practice — branch members map to the rep, which is set in step 5).
    for (int sweep = 0; sweep < 4; ++sweep) {
        bool changed = false;
        for (const auto& kv : bubbles.alt_map) {
            const uint32_t k = kv.first;
            const uint32_t t = kv.second;
            if (graph.is_phase_locked(k)) continue;
            const int8_t pt = graph.get_phase(t);
            if (pt != 0 && pt != 1) continue;
            graph.lock_phase(k, static_cast<int8_t>(1 - pt));
            changed = true;
        }
        if (!changed) break;
    }
}
