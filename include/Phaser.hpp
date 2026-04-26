#pragma once

#include <cstdint>
#include <memory>
#include <unordered_set>
#include <vector>

#include <thread_pool/thread_pool.hpp>

#include "BubbleDetector.hpp"   // BubbleResult
#include "TrioBinner.hpp"       // TrioScores

class Graph;
class ContactMatrix;

struct PhaserConfig {
    int      core_iterations = 200;
    int      sample_size     = 30;
    int      n_rounds        = 2;
    double   merge_threshold = 0.80;
    double   lock_ratio      = 0.80;
    uint32_t min_total_kmers = 1;
    uint32_t rng_seed        = 42;
};

size_t apply_trio_constraints(Graph& graph,
                               const std::unordered_set<uint32_t>& phasing_nodes,
                               const std::vector<TrioScores>&      scores,
                               const PhaserConfig&                 cfg = {});

void monte_carlo_phase(Graph&                                   graph,
                       const ContactMatrix&                     contacts,
                       const BubbleResult&                      bubbles,
                       const PhaserConfig&                      cfg,
                       std::shared_ptr<thread_pool::ThreadPool> pool);
