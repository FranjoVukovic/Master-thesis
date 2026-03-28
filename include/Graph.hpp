#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <utility>
#include <unordered_map>

struct Edge {
    uint32_t target_id;
    uint32_t overlap_length;
    bool target_rev;
};

// Hot metadata accessed by MCMC and phasing loops — kept small for cache efficiency.
struct NodeMetadata {
    int8_t   phase     = -1;   // -1 = unphased, 0 = haplotype 0, 1 = haplotype 1
    uint32_t length    = 0;
    uint32_t coverage  = 0;
    bool     is_reverse = false;
};

class Graph {
private:
    std::vector<uint64_t>      offsets;
    std::vector<Edge>          edges;
    std::vector<NodeMetadata>  node_info;   // hot: small fixed-size structs
    std::vector<std::string>   sequences_;  // cold: one entry per node, empty if GFA used '*'
    std::unordered_map<std::string, uint32_t> name_to_id;
    std::vector<std::string>   id_to_name;

public:
    Graph();
        
    void load_from_gfa(const std::string& filepath);

    // Optional: populate sequences_ from a companion FASTA file.
    // Call after load_from_gfa(). Only fills in nodes whose sequence is currently empty.
    // Needed when the GFA uses '*' for sequences (common with hifiasm / Shasta output).
    void load_fasta(const std::string& filepath);

    std::pair<std::vector<Edge>::const_iterator, std::vector<Edge>::const_iterator>
    get_neighbors(uint32_t node_id) const;

    uint32_t            get_id(const std::string& name) const;
    std::string         get_name(uint32_t id) const;
    const std::string&  get_sequence(uint32_t id) const;  // empty string if not available

    size_t get_num_nodes() const;
};