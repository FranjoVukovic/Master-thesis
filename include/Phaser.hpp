#pragma once

#include <cstdint>
#include <unordered_set>
#include <vector>

#include "TrioBinner.hpp"   // TrioScores

class Graph;
class ContactMatrix;

// ============================================================
// Trio-constrained initialization (data preparation stage).
//
// apply_trio_constraints():
//   For every node in phasing_nodes, inspects its TrioScores and — when one
//   parent's k-mer count dominates by >= lock_ratio of the total matched
//   k-mers — calls graph.lock_phase(id, 0 or 1).  Pat-dominant → phase 0,
//   mat-dominant → phase 1.  Ambiguous or zero-score nodes are left with
//   phase = -1 and is_phase_locked = false for a later phasing stage.
//   Returns the number of nodes that got locked.
// ============================================================

struct PhaserConfig {
    // Fraction of total matched (pat + mat) k-mers that one parent must
    // contribute to lock a node's phase.  0.8 = 80% as the brief specifies.
    double   lock_ratio = 0.8;
    // Minimum total matched k-mers before we consider a node lockable — below
    // this, any ratio is statistical noise.
    uint32_t min_total_kmers = 1;
};

// Applies trio binning constraints in-place on `graph`. Only nodes present
// in `phasing_nodes` are considered; returns how many got locked.
size_t apply_trio_constraints(Graph& graph,
                               const std::unordered_set<uint32_t>& phasing_nodes,
                               const std::vector<TrioScores>&      scores,
                               const PhaserConfig&                 cfg = {});
