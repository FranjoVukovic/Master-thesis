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

struct NodeMetadata {
    int8_t phase = -1;
    uint32_t length = 0;
    uint32_t coverage = 0;
    bool is_reverse = false;
    std::string sequence;  // raw nucleotide sequence from GFA S-line ('*' stored as empty)
};

class Graph {
private:
    std::vector<uint64_t> offsets;
    std::vector<Edge> edges;
    std::vector<NodeMetadata> node_info;
    std::unordered_map<std::string, uint32_t> name_to_id;
    std::vector<std::string> id_to_name;

public:
    Graph();
        
    void load_from_gfa(const std::string& filepath);

    std::pair<std::vector<Edge>::const_iterator, std::vector<Edge>::const_iterator> 
    get_neighbors(uint32_t node_id) const;

    uint32_t get_id(const std::string& name) const;
    std::string get_name(uint32_t id) const;
    const std::string& get_sequence(uint32_t id) const;

    size_t get_num_nodes() const;
};