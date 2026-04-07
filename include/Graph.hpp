#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <iosfwd>
#include <utility>
#include <unordered_map>

struct Edge {
    uint32_t target_id;
    uint32_t overlap_length;
    bool source_rev;
    bool target_rev;
};

struct NodeMetadata {
    int8_t   phase           = -1;    // -1 = unphased, 0 = haplotype 0, 1 = haplotype 1
    uint32_t length          = 0;
    uint32_t coverage        = 0;
    bool     is_reverse      = false;
    bool     is_phase_locked = false;
};

class Graph {
private:
    std::vector<uint64_t>      offsets;
    std::vector<Edge>          edges;
    std::vector<NodeMetadata>  node_info;  
    std::vector<std::string>   sequences_; 
    std::unordered_map<std::string, uint32_t> name_to_id;
    std::vector<std::string>   id_to_name;
    void load_from_stream(std::istream& in);

public:
    Graph();
        
    void load_from_gfa(const std::string& filepath);

    void load_from_gfa_string(const std::string& content);

    void load_fasta(const std::string& filepath);

    std::pair<std::vector<Edge>::const_iterator, std::vector<Edge>::const_iterator>
    get_neighbors(uint32_t node_id) const;

    uint32_t            get_id(const std::string& name) const;
    std::string         get_name(uint32_t id) const;
    const std::string&  get_sequence(uint32_t id) const;

    size_t get_num_nodes() const;

    int8_t get_phase       (uint32_t id) const { return node_info[id].phase; }
    bool   is_phase_locked (uint32_t id) const { return node_info[id].is_phase_locked; }
    void   set_phase       (uint32_t id, int8_t phase)   { node_info[id].phase = phase; }
    void   lock_phase      (uint32_t id, int8_t phase)   {
        node_info[id].phase           = phase;
        node_info[id].is_phase_locked = true;
    }
};