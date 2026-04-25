#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <thread_pool/thread_pool.hpp>

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

    void load_contacts(const std::string& filepath, const Graph& graph,
                       int num_threads = 1);

    void add_contact(uint32_t u, uint32_t v, uint32_t weight = 1);

    void build_csr(size_t num_nodes,
                   std::shared_ptr<thread_pool::ThreadPool> pool = nullptr);

    void filter_to_phasing_nodes(const std::unordered_set<uint32_t>& phasing_nodes,
                                 size_t num_nodes);

    uint32_t get_contact(uint32_t node_a, uint32_t node_b) const;

    size_t   row_begin(uint32_t u)     const { return row_offsets[u];     }
    size_t   row_end  (uint32_t u)     const { return row_offsets[u + 1]; }
    uint32_t col_at   (size_t   idx)   const { return col_indices[idx];   }
    uint32_t value_at (size_t   idx)   const { return values[idx];        }
};