#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <thread_pool/thread_pool.hpp>

class Graph;

struct TrioScores {
    uint32_t pat_count = 0;  // k-mers matching paternal database
    uint32_t mat_count = 0;  // k-mers matching maternal database
};

// Assigns each graph node a paternal/maternal k-mer score using yak databases.
// Both databases must have been built with the same k-mer size (default: k=31).
// Nodes whose sequence is empty (GFA '*') receive scores of 0.
class TrioBinner {
public:
    // Opens pat_yak and mat_yak databases. Throws std::runtime_error on failure.
    TrioBinner(const std::string& pat_yak, const std::string& mat_yak, int k = 31);
    ~TrioBinner();

    // Returns one TrioScores per node (indexed by node ID).
    std::vector<TrioScores> compute_scores(const Graph& graph,
        std::shared_ptr<thread_pool::ThreadPool> pool = nullptr) const;

private:
    void* pat_db_;  // yak_ch_t* — paternal k-mer database
    void* mat_db_;  // yak_ch_t* — maternal k-mer database
    int   k_;

    // Count how many k-mers in seq appear in the given yak database.
    uint32_t count_matching_kmers(const std::string& seq, void* db) const;
};
