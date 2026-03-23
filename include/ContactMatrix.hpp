#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <string>
#include <utility>

class Graph;

struct PairHash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1); 
    }
};

class ContactMatrix {
private:
    std::unordered_map<std::pair<uint32_t, uint32_t>, uint32_t, PairHash> contact_builder;

    std::vector<uint64_t> row_offsets;
    std::vector<uint32_t> col_indices;
    std::vector<uint32_t> values;

public:
    ContactMatrix();

    void load_contacts(const std::string& filepath, const Graph& graph);

    void build_csr(size_t num_nodes);

    uint32_t get_contact(uint32_t node_a, uint32_t node_b) const;
};