#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <thread_pool/thread_pool.hpp>

class Graph;

struct BubbleResult {
    std::unordered_set<uint32_t>            phasing_nodes;
    std::unordered_map<uint32_t, uint32_t>  alt_map;
};

class BubbleDetector {
public:
    BubbleResult get_alts_from_shasta_names(const Graph& graph) const;

    BubbleResult find_unlabeled_alts(const Graph& graph,
                                     int    k                 = 16,
                                     int    sketch_size       = 1000,
                                     double jaccard_threshold = 0.2,
                                     std::shared_ptr<thread_pool::ThreadPool> pool = nullptr) const;

    // Parse a Shasta Assembly-Phased.csv to build bubble pairs directly from
    // the assembler's authoritative phasing. Rows with Ploidy==2 are grouped
    // by (Bubble chain, Position in bubble chain, Component); each group
    // must contain exactly two haplotypes (0 and 1).
    BubbleResult from_shasta_csv(const Graph& graph,
                                 const std::string& csv_path) const;

private:
    struct Bubble {
        uint32_t              source;
        uint32_t              sink;
        std::vector<uint32_t> branch_a;
        std::vector<uint32_t> branch_b; 
    };

    
    std::vector<Bubble> find_superbubbles(const Graph& graph) const;


    std::string branch_sequence(const std::vector<uint32_t>& branch,
                                const Graph& graph) const;


    std::vector<uint64_t> compute_sketch(const std::string& seq,
                                         int k,
                                         int sketch_size) const;

    double estimate_jaccard(const std::vector<uint64_t>& a,
                            const std::vector<uint64_t>& b) const;


    enum class AssemblerStyle {
        Classic,       // Shasta classic:  "000001F"       / "000001F_003"
        PRStyle,       // Shasta PR:       "PR.000001.0"   / "PR.000001.1"
        DotSuffix,     // Generic dot:     "contig1.0"     / "contig1.1"
        HifiAsmStyle,  // Hifiasm phased:  "h1tg000001l"   / "h2tg000001l"   (TODO)
        VerkkoStyle,   // Verkko phased:   "haplotype1-…"  / "haplotype2-…"  (TODO)
    };
    AssemblerStyle detect_style(const Graph& graph) const;


    std::string base_name(const std::string& name, AssemblerStyle style) const;


    std::vector<std::vector<uint32_t>> build_directed_graph(const Graph& graph) const;


    std::vector<int32_t> compute_topological_ranks(const std::vector<std::vector<uint32_t>>& adj) const;
};
