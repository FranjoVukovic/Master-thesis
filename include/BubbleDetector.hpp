#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

class Graph;

// Nodes kept for phasing and their pairwise alt-allele relationships.
// - phasing_nodes: every node that is an interior branch node of a bubble.
//   Only these are kept in the contact matrix (non-bubble nodes are noise).
// - alt_map[u] = v means u and v are alternate alleles of the same locus.
//   The mapping is bidirectional (alt_map[v] = u as well).
struct BubbleResult {
    std::unordered_set<uint32_t>            phasing_nodes;
    std::unordered_map<uint32_t, uint32_t>  alt_map;
};

class BubbleDetector {
public:
    // Mode: shasta names.
    // Groups node names by their base name (auto-detects naming convention).
    // Groups of exactly 2 → bubble pair.  Groups of 1 or >2 → skipped.
    // Supports three Shasta conventions (auto-detected from the majority of names):
    //   Classic  : "000001F" / "000001F_003"   (strip trailing _NNN suffix)
    //   PR-style : "PR.000001.0" / "PR.000001.1" (strip trailing .N suffix)
    //   Dot-suffix: "contig1.0" / "contig1.1"   (strip trailing .N suffix)
    BubbleResult get_alts_from_shasta_names(const Graph& graph) const;

    // Mode: minhash.
    // Finds superbubbles topologically using chain-tracing (degree-2 nodes),
    // then confirms each bubble pair with MinHash k-mer similarity.
    // Bubbles whose branch sequences are unavailable (all '*') are skipped.
    BubbleResult find_unlabeled_alts(const Graph& graph,
                                     int    k                 = 16,
                                     int    sketch_size       = 1000,
                                     double jaccard_threshold = 0.2) const;

private:
    // A topological superbubble: two paths (chains) from source to sink.
    // branch_a / branch_b hold the interior nodes of each path in order,
    // excluding source and sink.
    struct Bubble {
        uint32_t              source;
        uint32_t              sink;
        std::vector<uint32_t> branch_a;
        std::vector<uint32_t> branch_b;
    };

    // Trace degree-2 chains to find all 2-branch superbubbles.
    std::vector<Bubble> find_superbubbles(const Graph& graph) const;

    // Concatenate sequences of all nodes in a branch (in order).
    // Returns empty string if any node has no sequence.
    std::string branch_sequence(const std::vector<uint32_t>& branch,
                                const Graph& graph) const;

    // Compute a MinHash sketch (sorted vector of the 'sketch_size' smallest hashes).
    // Uses canonical k-mers (min of forward and reverse-complement).
    std::vector<uint64_t> compute_sketch(const std::string& seq,
                                         int k,
                                         int sketch_size) const;

    // Estimate Jaccard similarity from two sorted MinHash sketches.
    double estimate_jaccard(const std::vector<uint64_t>& a,
                            const std::vector<uint64_t>& b) const;

    // Shasta name helpers.
    enum class ShastaStyle { Classic, PRStyle, DotSuffix };
    ShastaStyle detect_style(const Graph& graph) const;

    // Returns the base name for a node name under the detected style,
    // or empty string if the name does not match the pattern.
    std::string base_name(const std::string& name, ShastaStyle style) const;
};
