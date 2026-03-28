#include "BubbleDetector.hpp"
#include "Graph.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <queue>
#include <regex>
#include <stdexcept>
#include <unordered_map>

// ============================================================
// Superbubble finding via degree-2 chain tracing
// ============================================================
//
// In a bidirectional assembly graph every L-line creates u↔v and v↔u,
// so the graph is effectively undirected.
//
// Key observation:
//   - Interior nodes of a bubble branch connect ONLY to their two path
//     neighbours (source/sink or adjacent chain nodes), giving degree = 2.
//   - Source and sink nodes have degree ≥ 3 because they also connect to
//     the surrounding unique (non-bubble) regions.
//
// Algorithm:
//   1. Trace every connected chain of degree-2 nodes.
//      The two endpoints of a chain are the first nodes in each direction
//      with degree ≠ 2.
//   2. Group chains by their endpoint pair {ep0, ep1}.
//   3. If exactly 2 chains share the same endpoint pair → superbubble.

std::vector<BubbleDetector::Bubble>
BubbleDetector::find_superbubbles(const Graph& graph) const {
    const size_t n = graph.get_num_nodes();

    // Build a plain adjacency list (vector<vector>) for fast iteration.
    std::vector<std::vector<uint32_t>> adj(n);
    for (uint32_t u = 0; u < n; ++u) {
        auto [beg, end] = graph.get_neighbors(u);
        for (auto it = beg; it != end; ++it)
            adj[u].push_back(it->target_id);
    }

    // Trace chains of degree-2 nodes.
    // Map: normalised endpoint pair → list of interior-node chains.
    using EndpointPair = std::pair<uint32_t, uint32_t>;
    std::unordered_map<uint64_t, std::vector<std::vector<uint32_t>>> ep_chains;

    auto encode_ep = [](uint32_t a, uint32_t b) -> uint64_t {
        // Pack two uint32_t into one uint64_t key (a ≤ b).
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(a) << 32) | b;
    };

    std::vector<bool> visited(n, false);

    for (uint32_t seed = 0; seed < n; ++seed) {
        if (visited[seed] || adj[seed].size() != 2) continue;

        // --- trace left (direction 0) ---
        std::vector<uint32_t> left;
        {
            uint32_t prev = seed;
            uint32_t curr = adj[seed][0];
            while (!visited[curr] && adj[curr].size() == 2) {
                visited[curr] = true;
                left.push_back(curr);
                uint32_t next = (adj[curr][0] == prev) ? adj[curr][1] : adj[curr][0];
                prev = curr;
                curr = next;
            }
            // curr is now the left endpoint (degree ≠ 2 or already visited chain)
            left.push_back(curr);  // temporarily store endpoint at end
        }
        uint32_t ep0 = left.back();
        left.pop_back();  // remove endpoint, keep only interior

        // --- trace right (direction 1) ---
        std::vector<uint32_t> right;
        {
            uint32_t prev = seed;
            uint32_t curr = adj[seed][1];
            while (!visited[curr] && adj[curr].size() == 2) {
                visited[curr] = true;
                right.push_back(curr);
                uint32_t next = (adj[curr][0] == prev) ? adj[curr][1] : adj[curr][0];
                prev = curr;
                curr = next;
            }
            right.push_back(curr);
        }
        uint32_t ep1 = right.back();
        right.pop_back();

        visited[seed] = true;

        // Endpoints must be distinct (avoid self-loops collapsing to a single node).
        if (ep0 == ep1) continue;

        // Build the full interior chain: reversed(left) + seed + right.
        std::vector<uint32_t> chain;
        chain.reserve(left.size() + 1 + right.size());
        for (auto it = left.rbegin(); it != left.rend(); ++it)
            chain.push_back(*it);
        chain.push_back(seed);
        chain.insert(chain.end(), right.begin(), right.end());

        ep_chains[encode_ep(ep0, ep1)].push_back(std::move(chain));
    }

    // Collect bubbles: endpoint pairs with exactly 2 chains.
    std::vector<Bubble> bubbles;
    for (auto& [key, chains] : ep_chains) {
        if (chains.size() != 2) continue;  // only diploid (2-branch) bubbles

        uint32_t ep0 = static_cast<uint32_t>(key >> 32);
        uint32_t ep1 = static_cast<uint32_t>(key & 0xFFFFFFFF);

        bubbles.push_back({ep0, ep1,
                           std::move(chains[0]),
                           std::move(chains[1])});
    }

    return bubbles;
}

// ============================================================
// Branch sequence helpers
// ============================================================

std::string BubbleDetector::branch_sequence(const std::vector<uint32_t>& branch,
                                             const Graph& graph) const {
    std::string seq;
    for (uint32_t id : branch) {
        const std::string& s = graph.get_sequence(id);
        if (s.empty()) return {};  // node has no sequence
        seq += s;
    }
    return seq;
}

// ============================================================
// MinHash sketch
// ============================================================
//
// We use canonical k-mers (min of forward and rev-complement 2-bit encodings)
// hashed with a MurmurHash3-style 64-bit finaliser.
// The sketch is the 'sketch_size' smallest hash values (sorted).

static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

// 2-bit encoding: A=0, C=1, G=2, T=3.  Returns 4 for ambiguous bases.
static const uint8_t BASE2BIT[256] = {
    4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,
    4,4,4,4,4,4,4,4, 4,4,4,4,4,4,4,4,
    4,0,4,1,4,4,4,2, 4,4,4,4,4,4,4,4,   // A C G
    4,4,4,4,3,4,4,4, 4,4,4,4,4,4,4,4,   // T
    4,0,4,1,4,4,4,2, 4,4,4,4,4,4,4,4,   // a c g
    4,4,4,4,3,4,4,4, 4,4,4,4,4,4,4,4,   // t
};

std::vector<uint64_t> BubbleDetector::compute_sketch(const std::string& seq,
                                                      int k,
                                                      int sketch_size) const {
    if (static_cast<int>(seq.size()) < k) return {};

    const uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1) : UINT64_MAX;

    // We maintain a max-heap of size sketch_size holding the smallest hashes seen.
    std::vector<uint64_t> heap;
    heap.reserve(sketch_size + 1);

    uint64_t fwd = 0, rev = 0;
    int valid = 0;

    for (char c : seq) {
        uint8_t b = BASE2BIT[static_cast<unsigned char>(c)];
        if (b >= 4) { fwd = rev = 0; valid = 0; continue; }

        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 ^ b) << (2 * (k - 1)));
        rev &= mask;
        if (++valid < k) continue;

        uint64_t h = mix64(std::min(fwd, rev));  // canonical k-mer

        if (static_cast<int>(heap.size()) < sketch_size) {
            heap.push_back(h);
            if (static_cast<int>(heap.size()) == sketch_size)
                std::make_heap(heap.begin(), heap.end());  // max-heap
        } else if (h < heap.front()) {
            std::pop_heap(heap.begin(), heap.end());
            heap.back() = h;
            std::push_heap(heap.begin(), heap.end());
        }
    }

    std::sort(heap.begin(), heap.end());
    return heap;
}

double BubbleDetector::estimate_jaccard(const std::vector<uint64_t>& a,
                                         const std::vector<uint64_t>& b) const {
    if (a.empty() || b.empty()) return 0.0;

    // Count intersection of two sorted vectors.
    size_t shared = 0;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if      (a[i] < b[j]) ++i;
        else if (b[j] < a[i]) ++j;
        else                  { ++shared; ++i; ++j; }
    }

    // MinHash Jaccard estimate: shared / union-size (bounded by sketch).
    size_t union_size = a.size() + b.size() - shared;
    return union_size == 0 ? 0.0 : static_cast<double>(shared) / union_size;
}

// ============================================================
// find_unlabeled_alts  (topology + MinHash)
// ============================================================

BubbleResult BubbleDetector::find_unlabeled_alts(const Graph& graph,
                                                   int    k,
                                                   int    sketch_size,
                                                   double jaccard_threshold) const {
    BubbleResult result;

    for (const Bubble& b : find_superbubbles(graph)) {
        // Concatenate sequences for each branch.
        std::string seq_a = branch_sequence(b.branch_a, graph);
        std::string seq_b = branch_sequence(b.branch_b, graph);

        // If sequences are unavailable, skip MinHash but still record as phasing nodes
        // (the topology alone already confirms the bubble structure).
        bool seq_ok = !seq_a.empty() && !seq_b.empty();

        if (seq_ok) {
            auto sketch_a = compute_sketch(seq_a, k, sketch_size);
            auto sketch_b = compute_sketch(seq_b, k, sketch_size);
            double jac = estimate_jaccard(sketch_a, sketch_b);
            if (jac < jaccard_threshold) continue;  // sequences too dissimilar — skip
        }

        // Record phasing nodes and alt-map for both branches.
        uint32_t rep_a = b.branch_a.empty() ? b.source : b.branch_a[0];
        uint32_t rep_b = b.branch_b.empty() ? b.sink   : b.branch_b[0];

        for (uint32_t v : b.branch_a) {
            result.phasing_nodes.insert(v);
            result.alt_map[v] = rep_b;
        }
        for (uint32_t v : b.branch_b) {
            result.phasing_nodes.insert(v);
            result.alt_map[v] = rep_a;
        }
    }

    return result;
}

// ============================================================
// Shasta naming
// ============================================================

BubbleDetector::ShastaStyle
BubbleDetector::detect_style(const Graph& graph) const {
    // Classic  : name matches /^\d+[FR](_\d+)?$/   alt suffix is _NNN
    // PR-style : name matches /^PR\.\d+\.\d+$/      alt suffix is .\d
    // DotSuffix: name matches /^.+\.\d+$/            alt suffix is .\d
    //
    // Vote on a sample of up to 1000 node names.
    const std::regex re_classic  {R"(^\d+[FR](_\d+)?$)"};
    const std::regex re_pr       {R"(^PR\.\d+\.\d+$)"};
    const std::regex re_dot      {R"(^.+\.\d+$)"};

    int votes[3] = {0, 0, 0};
    size_t limit = std::min(graph.get_num_nodes(), size_t{1000});

    for (size_t i = 0; i < limit; ++i) {
        const std::string& name = graph.get_name(static_cast<uint32_t>(i));
        if (std::regex_match(name, re_classic)) ++votes[0];
        if (std::regex_match(name, re_pr))      ++votes[1];
        if (std::regex_match(name, re_dot))      ++votes[2];
    }

    int best = static_cast<int>(
        std::max_element(votes, votes + 3) - votes);

    switch (best) {
        case 0:  return ShastaStyle::Classic;
        case 1:  return ShastaStyle::PRStyle;
        default: return ShastaStyle::DotSuffix;
    }
}

std::string BubbleDetector::base_name(const std::string& name,
                                       ShastaStyle style) const {
    switch (style) {
        case ShastaStyle::Classic: {
            // Strip trailing _NNN  (e.g. "000001F_003" → "000001F")
            auto pos = name.rfind('_');
            if (pos != std::string::npos) {
                bool all_digits = true;
                for (size_t i = pos + 1; i < name.size(); ++i)
                    if (!std::isdigit(name[i])) { all_digits = false; break; }
                if (all_digits && pos + 1 < name.size())
                    return name.substr(0, pos);
            }
            return name;  // no suffix → this IS the primary
        }
        case ShastaStyle::PRStyle:
        case ShastaStyle::DotSuffix: {
            // Strip trailing .N  (e.g. "PR.000001.1" → "PR.000001", "contig1.1" → "contig1")
            auto pos = name.rfind('.');
            if (pos != std::string::npos && pos + 1 < name.size()) {
                bool all_digits = true;
                for (size_t i = pos + 1; i < name.size(); ++i)
                    if (!std::isdigit(name[i])) { all_digits = false; break; }
                if (all_digits)
                    return name.substr(0, pos);
            }
            return {};  // doesn't match — skip this node
        }
    }
    return {};
}

BubbleResult BubbleDetector::get_alts_from_shasta_names(const Graph& graph) const {
    ShastaStyle style = detect_style(graph);

    // Group node IDs by base name.
    std::unordered_map<std::string, std::vector<uint32_t>> groups;
    for (uint32_t id = 0; id < graph.get_num_nodes(); ++id) {
        std::string base = base_name(graph.get_name(id), style);
        if (!base.empty())
            groups[base].push_back(id);
    }

    BubbleResult result;
    for (auto& [base, ids] : groups) {
        if (ids.size() != 2) continue;  // only clean diploid pairs

        uint32_t a = ids[0], b = ids[1];
        result.phasing_nodes.insert(a);
        result.phasing_nodes.insert(b);
        result.alt_map[a] = b;
        result.alt_map[b] = a;
    }

    return result;
}
