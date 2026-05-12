#include "UnzippedGraph.hpp"

#include "Graph.hpp"

#include <stdexcept>

namespace {

bool is_homozygous(uint32_t id, const BubbleResult& bubbles) {
    return bubbles.phasing_nodes.find(id) == bubbles.phasing_nodes.end();
}

const char* phase_suffix(int8_t phase) {
    if (phase == 0) return "_ph0";
    if (phase == 1) return "_ph1";
    return "";
}

struct EdgeInfo {
    bool     found          = false;
    bool     source_rev     = false;
    bool     target_rev     = false;
    uint32_t overlap_length = 0;
};

EdgeInfo find_forward_edge(const Graph& g, uint32_t u, uint32_t v) {
    auto range = g.get_neighbors(u);
    for (auto it = range.first; it != range.second; ++it) {
        if (it->source_rev) continue;
        if (it->target_id == v) return {true, it->source_rev, it->target_rev, it->overlap_length};
    }
    // Fallback: any edge.
    for (auto it = range.first; it != range.second; ++it) {
        if (it->target_id == v) return {true, it->source_rev, it->target_rev, it->overlap_length};
    }
    return {};
}

void emit_chain(const Graph& g,
                const BubbleResult& bubbles,
                const Chain& chain,
                std::vector<UzSegment>& segs,
                std::vector<UzLink>& links) {
    std::vector<uint32_t> idxs;
    idxs.reserve(chain.nodes.size());

    for (uint32_t nid : chain.nodes) {
        UzSegment s;
        s.src_id   = nid;
        s.phase    = chain.phase;
        const uint32_t md_len = g.get_length(nid);
        s.length   = md_len ? md_len : static_cast<uint32_t>(g.get_sequence(nid).size());
        s.coverage = g.get_coverage(nid);

        std::string base = g.get_name(nid);
        if (is_homozygous(nid, bubbles) && chain.phase != -1) {
            s.name = base + phase_suffix(chain.phase);
        } else {
            s.name = base;
        }
        idxs.push_back(static_cast<uint32_t>(segs.size()));
        segs.push_back(std::move(s));
    }

    for (size_t i = 0; i + 1 < chain.nodes.size(); ++i) {
        EdgeInfo e = find_forward_edge(g, chain.nodes[i], chain.nodes[i + 1]);
        UzLink l;
        l.source_idx     = idxs[i];
        l.target_idx     = idxs[i + 1];
        l.source_rev     = e.source_rev;
        l.target_rev     = e.target_rev;
        l.overlap_length = e.overlap_length;
        links.push_back(l);
    }
}

} // namespace

UnzippedGraph UnzippedGraph::build(const Graph& graph,
                                   const BubbleResult& bubbles,
                                   const ChainResult& chains) {
    UnzippedGraph uz;
    for (const auto& c : chains.phase_0)  emit_chain(graph, bubbles, c, uz.segs_, uz.links_);
    for (const auto& c : chains.phase_1)  emit_chain(graph, bubbles, c, uz.segs_, uz.links_);
    for (const auto& c : chains.unphased) emit_chain(graph, bubbles, c, uz.segs_, uz.links_);
    return uz;
}
