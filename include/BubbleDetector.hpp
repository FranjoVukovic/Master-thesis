#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

class Graph;

struct BubbleResult {
    std::unordered_set<uint32_t>            phasing_nodes;
    std::unordered_map<uint32_t, uint32_t>  alt_map;
};

class BubbleDetector {
public:
    // Currently implemented conventions (auto-detected from the majority of names):
    //   Classic   : "000001F" / "000001F_003"     (strip trailing _NNN suffix)
    //   PR-style  : "PR.000001.0" / "PR.000001.1" (strip trailing .N suffix)
    //   DotSuffix : "contig1.0" / "contig1.1"     (strip trailing .N suffix)
    //   HifiAsm   : "h1tg000001l" / "h2tg000001l" — detected but NOT YET implemented
    //   Verkko    : "haplotype1-0000001" / "haplotype2-0000001" — detected but NOT YET implemented
    BubbleResult get_alts_from_shasta_names(const Graph& graph) const;

    BubbleResult find_unlabeled_alts(const Graph& graph,
                                     int    k                 = 16,
                                     int    sketch_size       = 1000,
                                     double jaccard_threshold = 0.2) const;

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
