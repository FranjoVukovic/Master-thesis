#include "Graph.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

struct RawLink {
    uint32_t source_id;
    uint32_t target_id;
    uint32_t overlap_length;
    bool source_rev;
    bool target_rev;
};

uint32_t parse_cigar_length(const std::string& cigar) {
    if (cigar == "*" || cigar.empty()) return 0;
    
    uint32_t total_length = 0;
    uint32_t current_num = 0;
    
    for (char c : cigar) {
        if (std::isdigit(c)) {
            current_num = current_num * 10 + (c - '0');
        } else {
            total_length += current_num;
            current_num = 0;
        }
    }
    return total_length;
}

Graph::Graph() = default;

void Graph::load_from_gfa(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open())
        throw std::runtime_error("Cannot open GFA file: " + filepath);
    load_from_stream(file);
}

void Graph::load_from_gfa_string(const std::string& content) {
    std::istringstream iss(content);
    load_from_stream(iss);
}

void Graph::load_from_stream(std::istream& in) {
    std::string line;
    std::vector<RawLink> temp_links;
    uint32_t current_id = 0;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == 'S') {
            std::istringstream iss(line);
            std::string type, name, sequence, token;
            iss >> type >> name >> sequence;

            std::string seq_stored = (sequence == "*") ? "" : sequence;
            uint32_t node_length = seq_stored.empty() ? 0 : static_cast<uint32_t>(seq_stored.size());
            uint32_t node_coverage = 0;

            while (iss >> token) {
                if (token.rfind("LN:i:", 0) == 0) {
                    node_length = std::stoul(token.substr(5));
                } else if (token.rfind("RC:i:", 0) == 0) {
                    node_coverage = std::stoul(token.substr(5));
                }
            }

            if (name_to_id.find(name) == name_to_id.end()) {
                name_to_id[name] = current_id++;
                id_to_name.push_back(name);

                NodeMetadata meta;
                meta.length   = node_length;
                meta.coverage = node_coverage;
                node_info.push_back(meta);
                sequences_.push_back(std::move(seq_stored));
            } else {
                uint32_t id = name_to_id[name];
                node_info[id].length   = node_length;
                node_info[id].coverage = node_coverage;
                sequences_[id] = std::move(seq_stored);
            }
        } 
        else if (line[0] == 'L') {
            std::istringstream iss(line);
            std::string type, src_name, src_orient, tgt_name, tgt_orient, overlap;
            iss >> type >> src_name >> src_orient >> tgt_name >> tgt_orient >> overlap;

            if (name_to_id.find(src_name) == name_to_id.end()) {
                name_to_id[src_name] = current_id++;
                id_to_name.push_back(src_name);
                node_info.push_back(NodeMetadata());
                sequences_.emplace_back();
            }
            if (name_to_id.find(tgt_name) == name_to_id.end()) {
                name_to_id[tgt_name] = current_id++;
                id_to_name.push_back(tgt_name);
                node_info.push_back(NodeMetadata());
                sequences_.emplace_back();
            }

            uint32_t u = name_to_id[src_name];
            uint32_t v = name_to_id[tgt_name];
            bool u_rev = (src_orient == "-");
            bool v_rev = (tgt_orient == "-");

            uint32_t overlap_length = 0;
            if (!overlap.empty() && overlap != "*" ) {
                overlap_length = parse_cigar_length(overlap);
            }

            temp_links.push_back({u, v, overlap_length, u_rev, v_rev});    
            temp_links.push_back({v, u, overlap_length, !v_rev, u_rev});
        }
    }

    uint32_t num_nodes = current_id;
    offsets.assign(num_nodes + 1, 0);

    for (const auto& link : temp_links) {
        offsets[link.source_id + 1]++;
    }

    for (size_t i = 1; i <= num_nodes; ++i) {
        offsets[i] += offsets[i - 1];
    }

    edges.resize(temp_links.size());
    
    std::vector<uint64_t> current_offsets = offsets; 

    for (const auto& link : temp_links) {
        uint64_t pos = current_offsets[link.source_id]++;
        edges[pos] = {link.target_id, link.overlap_length, link.source_rev, link.target_rev};
    }
}

std::pair<std::vector<Edge>::const_iterator, std::vector<Edge>::const_iterator> 
Graph::get_neighbors(uint32_t node_id) const {
    auto start = edges.begin() + offsets[node_id];
    auto end = edges.begin() + offsets[node_id + 1];
    return {start, end};
}

uint32_t Graph::get_id(const std::string& name) const {
    return name_to_id.at(name);
}

std::string Graph::get_name(uint32_t id) const {
    return id_to_name[id];
}

const std::string& Graph::get_sequence(uint32_t id) const {
    return sequences_[id];
}

void Graph::load_fasta(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open())
        throw std::runtime_error("Cannot open FASTA file: " + filepath);

    std::string line, current_name, current_seq;

    auto flush = [&]() {
        if (current_name.empty()) return;
        auto it = name_to_id.find(current_name);
        if (it != name_to_id.end()) {
            uint32_t id = it->second;
            if (sequences_[id].empty())
                sequences_[id] = std::move(current_seq);
        }
        current_name.clear();
        current_seq.clear();
    };

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            flush();
            size_t space = line.find_first_of(" \t", 1);
            current_name = line.substr(1, space - 1);
        } else {
            current_seq += line;
        }
    }
    flush();
}

size_t Graph::get_num_nodes() const {
    return id_to_name.size();
}