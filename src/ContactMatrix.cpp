#include "ContactMatrix.hpp"
#include "Graph.hpp"
#include "BamParser.hpp" 
#include <algorithm>
#include <stdexcept>
#include <iostream>

ContactMatrix::ContactMatrix() = default;

void ContactMatrix::load_contacts(const std::string& filepath, const Graph& graph) {
    BamParser parser(filepath);
    std::string name_a, name_b;

    while (parser.get_next_contact(name_a, name_b)) {
        try {
            uint32_t u = graph.get_id(name_a);
            uint32_t v = graph.get_id(name_b);
            
            if (u > v) {
                std::swap(u, v);
            }
            
            contact_builder[{u, v}] += 1;
            
        } catch (const std::out_of_range&) {
            continue;
        }
    }
}

void ContactMatrix::add_contact(uint32_t u, uint32_t v, uint32_t weight) {
    if (u > v) std::swap(u, v);
    contact_builder[{u, v}] += weight;
}

void ContactMatrix::build_csr(size_t num_nodes) {
    row_offsets.assign(num_nodes + 1, 0);

    for (const auto& kv : contact_builder) {
        uint32_t u = kv.first.first;
        uint32_t v = kv.first.second;
        
        row_offsets[u + 1]++;
        if (u != v) {
            row_offsets[v + 1]++;
        }
    }

    for (size_t i = 0; i < num_nodes; ++i) {
        row_offsets[i + 1] += row_offsets[i];
    }

    size_t total_elements = row_offsets.back();
    col_indices.resize(total_elements);
    values.resize(total_elements);

    std::vector<uint64_t> current_offsets = row_offsets;

    for (const auto& kv : contact_builder) {
        uint32_t u = kv.first.first;
        uint32_t v = kv.first.second;
        uint32_t weight = kv.second;

        uint64_t pos_u = current_offsets[u]++;
        col_indices[pos_u] = v;
        values[pos_u] = weight;

        if (u != v) {
            uint64_t pos_v = current_offsets[v]++;
            col_indices[pos_v] = u;
            values[pos_v] = weight;
        }
    }

    contact_builder.clear();

    for (size_t i = 0; i < num_nodes; ++i) {
        uint64_t row_start = row_offsets[i];
        uint64_t row_end   = row_offsets[i + 1];
        if (row_end - row_start <= 1) continue;

        std::vector<std::pair<uint32_t,uint32_t>> row_data;
        row_data.reserve(row_end - row_start);
        for (uint64_t j = row_start; j < row_end; ++j)
            row_data.push_back({col_indices[j], values[j]});
        std::sort(row_data.begin(), row_data.end());
        for (size_t j = 0; j < row_data.size(); ++j) {
            col_indices[row_start + j] = row_data[j].first;
            values[row_start + j]      = row_data[j].second;
        }
    }
}

void ContactMatrix::filter_to_phasing_nodes(const std::unordered_set<uint32_t>& phasing_nodes,
                                             size_t num_nodes) {
    std::vector<uint64_t> new_offsets(num_nodes + 1, 0);

    for (size_t u = 0; u < num_nodes; ++u) {
        if (!phasing_nodes.count(static_cast<uint32_t>(u))) continue;
        uint64_t start = row_offsets[u];
        uint64_t end   = row_offsets[u + 1];
        for (uint64_t j = start; j < end; ++j) {
            if (phasing_nodes.count(col_indices[j]))
                new_offsets[u + 1]++;
        }
    }

    for (size_t i = 0; i < num_nodes; ++i)
        new_offsets[i + 1] += new_offsets[i];

    std::vector<uint32_t> new_cols(new_offsets.back());
    std::vector<uint32_t> new_vals(new_offsets.back());
    std::vector<uint64_t> cur = new_offsets;

    for (size_t u = 0; u < num_nodes; ++u) {
        if (!phasing_nodes.count(static_cast<uint32_t>(u))) continue;
        uint64_t start = row_offsets[u];
        uint64_t end   = row_offsets[u + 1];
        for (uint64_t j = start; j < end; ++j) {
            if (phasing_nodes.count(col_indices[j])) {
                uint64_t pos = cur[u]++;
                new_cols[pos] = col_indices[j];
                new_vals[pos] = values[j];
            }
        }
    }

    row_offsets = std::move(new_offsets);
    col_indices = std::move(new_cols);
    values      = std::move(new_vals);
}

uint32_t ContactMatrix::get_contact(uint32_t node_a, uint32_t node_b) const {
    uint64_t start = row_offsets[node_a];
    uint64_t end   = row_offsets[node_a + 1];

    auto it = std::lower_bound(col_indices.begin() + start,
                               col_indices.begin() + end,
                               node_b);
    if (it != col_indices.begin() + end && *it == node_b)
        return values[it - col_indices.begin()];

    return 0;
}