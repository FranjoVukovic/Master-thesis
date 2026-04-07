#include "Phaser.hpp"

#include "Graph.hpp"

#include <cstdint>
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
            graph.lock_phase(id, 0);  // paternal → haplotype 0
            ++locked;
        } else if (mat_frac >= cfg.lock_ratio) {
            graph.lock_phase(id, 1);  // maternal → haplotype 1
            ++locked;
        }
        // Ambiguous or low-total nodes stay phase = -1, is_phase_locked = false.
    }
    return locked;
}

