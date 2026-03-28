#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <string>
#include <utility>

class Graph;

struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1,T2>& p) const {
        std::size_t h = std::hash<T1>{}(p.first);
        h ^= std::hash<T2>{}(p.second) + 0x9e3779b9u + (h << 6) + (h >> 2);
        return h;
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

    // Remove all contacts involving nodes NOT in phasing_nodes.
    // Call after build_csr(). Rebuilds the CSR in-place.
    void filter_to_bubbles(const std::unordered_set<uint32_t>& phasing_nodes,
                           size_t num_nodes);

    uint32_t get_contact(uint32_t node_a, uint32_t node_b) const;
};