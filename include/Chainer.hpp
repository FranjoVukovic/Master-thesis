#pragma once

#include <cstdint>
#include <vector>

#include "BubbleDetector.hpp"

class Graph;
class ContactMatrix;

struct Chain {
    std::vector<uint32_t> nodes;
    std::vector<bool>     reverse;
    int8_t                phase = -1;
};

struct ChainResult {
    std::vector<Chain> phase_0;
    std::vector<Chain> phase_1;
    std::vector<Chain> unphased;
};

class Chainer {
public:
    ChainResult generate_chain_paths(const Graph& graph,
                                     const BubbleResult& bubbles,
                                     const ContactMatrix* contacts = nullptr) const;
};
