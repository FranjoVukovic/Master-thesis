#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "BubbleDetector.hpp"
#include "Chainer.hpp"

class Graph;

struct UzSegment {
    std::string name;
    uint32_t    length   = 0;
    uint32_t    coverage = 0;
    int8_t      phase    = -1;
    uint32_t    src_id   = 0;
};

struct UzLink {
    uint32_t source_idx     = 0;
    uint32_t target_idx     = 0;
    bool     source_rev     = false;
    bool     target_rev     = false;
    uint32_t overlap_length = 0;
};

class UnzippedGraph {
public:
    static UnzippedGraph build(const Graph& graph,
                               const BubbleResult& bubbles,
                               const ChainResult& chains);

    const std::vector<UzSegment>& segments() const { return segs_; }
    const std::vector<UzLink>&    links()    const { return links_; }

private:
    std::vector<UzSegment> segs_;
    std::vector<UzLink>    links_;
};
