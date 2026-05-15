#include "BubbleDetector.hpp"
#include "Graph.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <future>
#include <numeric>
#include <queue>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

std::vector<BubbleDetector::Bubble>
BubbleDetector::find_superbubbles(const Graph& graph) const {
    auto adj = build_directed_graph(graph);
    auto rank = compute_topological_ranks(adj);

    int32_t num_valid = 0;
    for (int32_t r : rank) {
        if (r >= num_valid) num_valid = r + 1;
    }

    if (num_valid < 4) return {};

    std::vector<uint32_t> rank_to_split(num_valid);
    for (uint32_t split_node = 0; split_node < rank.size(); ++split_node) {
        if (rank[split_node] != -1) {
            rank_to_split[rank[split_node]] = split_node;
        }
    }

    std::vector<int32_t> out_cnt(num_valid, 0), in_cnt(num_valid, 0);
    std::vector<int32_t> out_child(num_valid, -1); 
    std::vector<int32_t> in_parent(num_valid,  num_valid);  

    for (int32_t u_rank = 0; u_rank < num_valid; ++u_rank) {
        uint32_t u_split = rank_to_split[u_rank];
        for (uint32_t v_split : adj[u_split]) {
            int32_t v_rank = rank[v_split];
            if (v_rank == -1) continue; 

            if (u_rank < v_rank) {
                ++out_cnt[u_rank]; ++in_cnt[v_rank];
                if (v_rank > out_child[u_rank]) out_child[u_rank] = v_rank;
                if (u_rank < in_parent[v_rank]) in_parent[v_rank] = u_rank;
            }
        }
    }

    std::vector<uint32_t> out_off(num_valid + 1, 0);
    for (int32_t i = 0; i < num_valid; ++i) out_off[i + 1] = out_off[i] + out_cnt[i];
    std::vector<uint32_t> out_nbr(out_off[num_valid]);
    {
        std::vector<uint32_t> op(num_valid, 0);
        for (int32_t u_rank = 0; u_rank < num_valid; ++u_rank) {
            uint32_t u_split = rank_to_split[u_rank];
            for (uint32_t v_split : adj[u_split]) {
                int32_t v_rank = rank[v_split];
                if (v_rank != -1 && u_rank < v_rank) {
                    out_nbr[out_off[u_rank] + op[u_rank]++] = static_cast<uint32_t>(v_rank);
                }
            }
        }
    }

    int LOG = 1;
    while ((1 << LOG) <= num_valid) ++LOG;

    std::vector<std::vector<int32_t>> sp_max(LOG, std::vector<int32_t>(num_valid, -1));
    std::vector<std::vector<int32_t>> sp_min(LOG, std::vector<int32_t>(num_valid,  num_valid));
    for (int32_t i = 0; i < num_valid; ++i) { sp_max[0][i] = out_child[i]; sp_min[0][i] = in_parent[i]; }
    for (int k = 1; k < LOG; ++k) {
        const int32_t half = 1 << (k - 1);
        for (int32_t i = 0; i + (1 << k) <= num_valid; ++i) {
            sp_max[k][i] = std::max(sp_max[k-1][i], sp_max[k-1][i + half]);
            sp_min[k][i] = std::min(sp_min[k-1][i], sp_min[k-1][i + half]);
        }
    }
    std::vector<int> lg(num_valid + 1, 0);
    for (int i = 2; i <= num_valid; ++i) lg[i] = lg[i / 2] + 1;

    auto rmax = [&](int32_t l, int32_t r) -> int32_t {
        if (l > r) return -1;
        int k = lg[r - l + 1];
        return std::max(sp_max[k][l], sp_max[k][r - (1 << k) + 1]);
    };
    auto rmin = [&](int32_t l, int32_t r) -> int32_t {
        if (l > r) return num_valid;
        int k = lg[r - l + 1];
        return std::min(sp_min[k][l], sp_min[k][r - (1 << k) + 1]);
    };

    std::vector<bool> is_exit(num_valid, false), is_entrance(num_valid, false);

    std::vector<std::vector<uint32_t>> in_nbrs(num_valid);
    for (int32_t u = 0; u < num_valid; ++u) {
        for (uint32_t j = out_off[u]; j < out_off[u + 1]; ++j) {
            in_nbrs[out_nbr[j]].push_back(u);
        }
    }

    for (int32_t v_rank = 0; v_rank < num_valid; ++v_rank) {
        for (uint32_t u : in_nbrs[v_rank]) {
            if (out_cnt[u] == 1) { 
                is_exit[v_rank] = true; 
                break; 
            }
        }

        for (uint32_t j = out_off[v_rank]; j < out_off[v_rank + 1]; ++j) {
            if (in_cnt[static_cast<int32_t>(out_nbr[j])] == 1) { 
                is_entrance[v_rank] = true; 
                break; 
            }
        }
    }

    struct Cand { int32_t node; bool is_ent; }; 
    std::vector<Cand> cands;
    cands.reserve(2 * num_valid);

    std::vector<int32_t> prev_ent(num_valid, -1);
    std::vector<int32_t> alt_ent(num_valid,  -1);
    std::vector<int32_t> ent_pos(num_valid,  -1);

    int32_t last_ent = -1;
    for (int32_t v = 0; v < num_valid; ++v) {
        prev_ent[v] = last_ent;
        if (is_exit[v])     cands.push_back({v, false});
        if (is_entrance[v]) { ent_pos[v] = static_cast<int32_t>(cands.size()); cands.push_back({v, true}); last_ent = v; }
    }

    auto validate = [&](int32_t s, int32_t t) -> int32_t {
        if (s < 0 || t <= s + 1) return -1;
        if (rmax(s, t - 1) != t) return -1;
        int32_t op = rmin(s + 1, t);
        if (op == s)   return s;
        if (op >= 0 && is_entrance[op]) return op;
        if (op >= 0) return prev_ent[op];
        return -1;
    };

    auto bfs_branch = [&](uint32_t head, uint32_t s_id, uint32_t t_id)
            -> std::vector<uint32_t> {
        const uint32_t span = t_id - s_id - 1;
        std::vector<uint32_t> res;
        if (span == 0) return res;
        std::vector<uint8_t> vis(span, 0);
        res.push_back(head); vis[head - s_id - 1] = 1;
        for (size_t qi = 0; qi < res.size(); ++qi) {
            for (uint32_t j = out_off[res[qi]]; j < out_off[res[qi] + 1]; ++j) {
                uint32_t w = out_nbr[j];
                if (w > s_id && w < t_id && !vis[w - s_id - 1]) {
                    vis[w - s_id - 1] = 1; res.push_back(w);
                }
            }
        }
        return res;
    };

    std::vector<Bubble> bubbles;

    auto rank_to_phys = [&](const std::vector<uint32_t>& ranks) {
        std::vector<uint32_t> phys;
        phys.reserve(ranks.size());
        for (uint32_t r : ranks) {
            phys.push_back(rank_to_split[r] / 2); 
        }
        std::sort(phys.begin(), phys.end());
        phys.erase(std::unique(phys.begin(), phys.end()), phys.end());
        return phys;
    };

    auto emit = [&](int32_t s, int32_t t) {
        uint32_t su = static_cast<uint32_t>(s), tu = static_cast<uint32_t>(t);
        std::vector<uint32_t> heads;
        for (uint32_t j = out_off[su]; j < out_off[su + 1]; ++j) {
            uint32_t w = out_nbr[j];
            if (w > su && w < tu) heads.push_back(w);
        }
        if (heads.size() < 2) return;
        std::vector<uint32_t> ba, bb;
        if (heads.size() == 2) {
            ba = bfs_branch(heads[0], su, tu);
            bb = bfs_branch(heads[1], su, tu);
        } else {
            ba.reserve(tu - su - 1);
            for (uint32_t id = su + 1; id < tu; ++id) ba.push_back(id);
        }
        
        uint32_t phys_s = rank_to_split[su] / 2;
        uint32_t phys_t = rank_to_split[tu] / 2;

        bubbles.push_back({phys_s, phys_t, rank_to_phys(ba), rank_to_phys(bb)});
    };

    int32_t tail = static_cast<int32_t>(cands.size()) - 1;

    std::function<void(int32_t)> rsb = [&](int32_t start_node) {
        if (tail < 0) return;
        
        int32_t exit_node = cands[tail].node;

        if (start_node < 0 || exit_node < 0 || start_node >= exit_node) {
            --tail; 
            return;
        }

        int32_t s     = prev_ent[exit_node];
        int32_t valid = -1;
        while (s >= start_node && s >= 0) {
            valid = validate(s, exit_node);
            if (valid == s || valid == alt_ent[s] || valid == -1) break;
            alt_ent[s] = valid;
            s = valid;
        }

        --tail;   

        if (s >= start_node && s >= 0 && valid == s) {
            emit(s, exit_node);

            const int32_t sp = ent_pos[s]; 

            while (tail > sp) {
                if (!cands[tail].is_ent) {
                    const int32_t nxt = (sp + 1 <= tail) ? cands[sp + 1].node : -1;
                    rsb(nxt); 
                } else {
                    --tail;   
                }
            }
        }
    };

    while (tail >= 0) {
        if (cands[tail].is_ent) {
            --tail; 
        } else {
            rsb(cands[0].node);
        }
    }

    return bubbles;
}

std::string BubbleDetector::branch_sequence(const std::vector<uint32_t>& branch,
                                             const Graph& graph) const {
    std::string seq;
    for (uint32_t id : branch) {
        const std::string& s = graph.get_sequence(id);
        if (s.empty()) continue;  
        seq += s;
    }
    return seq;  
}

std::vector<std::vector<uint32_t>> BubbleDetector::build_directed_graph(const Graph& graph) const {
    const uint32_t n = static_cast<uint32_t>(graph.get_num_nodes());
    
    std::vector<std::vector<uint32_t>> adj(2 * n);

    for (uint32_t u = 0; u < n; ++u) {
        auto [beg, end] = graph.get_neighbors(u);
        for (auto it = beg; it != end; ++it) {
            uint32_t v = it->target_id;
            
            uint32_t u_split = it->source_rev ? (2 * u + 1) : (2 * u);
            uint32_t v_split = it->target_rev ? (2 * v + 1) : (2 * v);
            
            adj[u_split].push_back(v_split);
        }
    }
    
    return adj;
}

std::vector<int32_t> BubbleDetector::compute_topological_ranks(const std::vector<std::vector<uint32_t>>& adj) const {
    const uint32_t n = static_cast<uint32_t>(adj.size());
    enum Color : uint8_t { WHITE = 0, GRAY = 1, BLACK = 2 };
    std::vector<uint8_t>  color(n, WHITE);
    std::vector<uint32_t> post_order;
    post_order.reserve(n);
    std::vector<std::pair<uint32_t, size_t>> stack;
    stack.reserve(n);

    for (uint32_t root = 0; root < n; ++root) {
        if (color[root] != WHITE) continue;
        color[root] = GRAY;
        stack.push_back({root, 0});
        while (!stack.empty()) {
            auto& frame = stack.back();
            const uint32_t u = frame.first;
            if (frame.second < adj[u].size()) {
                const uint32_t v = adj[u][frame.second++];
                if (color[v] == WHITE) {
                    color[v] = GRAY;
                    stack.push_back({v, 0});
                }
            } else {
                color[u] = BLACK;
                post_order.push_back(u);
                stack.pop_back();
            }
        }
    }

    std::vector<int32_t> rank(n, -1);
    int32_t current_rank = 0;
    for (auto it = post_order.rbegin(); it != post_order.rend(); ++it) {
        rank[*it] = current_rank++;
    }
    return rank;
}

static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

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

        uint64_t h = mix64(std::min(fwd, rev));  

        if (static_cast<int>(heap.size()) < sketch_size) {
            heap.push_back(h);
            if (static_cast<int>(heap.size()) == sketch_size)
                std::make_heap(heap.begin(), heap.end());  
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

    size_t shared = 0;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if      (a[i] < b[j]) ++i;
        else if (b[j] < a[i]) ++j;
        else                  { ++shared; ++i; ++j; }
    }

    size_t union_size = a.size() + b.size() - shared;
    return union_size == 0 ? 0.0 : static_cast<double>(shared) / union_size;
}


BubbleResult BubbleDetector::find_unlabeled_alts(const Graph& graph,
                                                   int    k,
                                                   int    sketch_size,
                                                   double jaccard_threshold,
                                                   std::shared_ptr<thread_pool::ThreadPool> pool) const {
    auto bubbles = find_superbubbles(graph);
    BubbleResult result;

    // Handle single-branch bubbles (no sketch needed)
    for (const Bubble& b : bubbles) {
        if (b.branch_b.empty()) {
            for (uint32_t v : b.branch_a) result.phasing_nodes.insert(v);
        }
    }

    // Collect two-branch bubbles for sketch computation
    std::vector<size_t> two_branch_indices;
    for (size_t i = 0; i < bubbles.size(); ++i) {
        if (!bubbles[i].branch_b.empty()) two_branch_indices.push_back(i);
    }

    if (!pool || two_branch_indices.empty()) {
        for (size_t idx : two_branch_indices) {
            const Bubble& b = bubbles[idx];
            std::string seq_a = branch_sequence(b.branch_a, graph);
            std::string seq_b = branch_sequence(b.branch_b, graph);
            if (!seq_a.empty() && !seq_b.empty()) {
                double jac = estimate_jaccard(compute_sketch(seq_a, k, sketch_size),
                                              compute_sketch(seq_b, k, sketch_size));
                if (jac < jaccard_threshold) continue;
            }
            uint32_t rep_a = b.branch_a[0], rep_b = b.branch_b[0];
            for (uint32_t v : b.branch_a) { result.phasing_nodes.insert(v); result.alt_map[v] = rep_b; }
            for (uint32_t v : b.branch_b) { result.phasing_nodes.insert(v); result.alt_map[v] = rep_a; }
        }
    } else {
        const size_t num_threads = pool->num_threads();
        const size_t chunk_size = std::max<size_t>(1, (two_branch_indices.size() + num_threads - 1) / num_threads);
        std::vector<std::future<BubbleResult>> futures;

        for (size_t start = 0; start < two_branch_indices.size(); start += chunk_size) {
            size_t end = std::min(start + chunk_size, two_branch_indices.size());
            futures.emplace_back(pool->Submit(
                [this, &graph, &bubbles, &two_branch_indices, k, sketch_size, jaccard_threshold]
                (size_t begin, size_t end) {
                    BubbleResult local;
                    for (size_t ci = begin; ci < end; ++ci) {
                        const Bubble& b = bubbles[two_branch_indices[ci]];
                        std::string seq_a = branch_sequence(b.branch_a, graph);
                        std::string seq_b = branch_sequence(b.branch_b, graph);
                        if (!seq_a.empty() && !seq_b.empty()) {
                            double jac = estimate_jaccard(compute_sketch(seq_a, k, sketch_size),
                                                          compute_sketch(seq_b, k, sketch_size));
                            if (jac < jaccard_threshold) continue;
                        }
                        uint32_t rep_a = b.branch_a[0], rep_b = b.branch_b[0];
                        for (uint32_t v : b.branch_a) { local.phasing_nodes.insert(v); local.alt_map[v] = rep_b; }
                        for (uint32_t v : b.branch_b) { local.phasing_nodes.insert(v); local.alt_map[v] = rep_a; }
                    }
                    return local;
                },
                start, end));
        }

        for (auto& f : futures) {
            BubbleResult local = f.get();
            result.phasing_nodes.insert(local.phasing_nodes.begin(), local.phasing_nodes.end());
            result.alt_map.insert(local.alt_map.begin(), local.alt_map.end());
        }
    }

    // Sibling-MinHash pass: find alts among non-phased nodes using shared parents/successors.
    const uint32_t n_phys = static_cast<uint32_t>(graph.get_num_nodes());
    std::vector<std::vector<uint32_t>> out_phys(n_phys), in_phys(n_phys);
    for (uint32_t u = 0; u < n_phys; ++u) {
        auto [beg, end] = graph.get_neighbors(u);
        for (auto it = beg; it != end; ++it) {
            if (!it->source_rev && !it->target_rev) {
                out_phys[u].push_back(it->target_id);
                in_phys[it->target_id].push_back(u);
            }
        }
    }

    std::unordered_map<uint32_t, std::vector<uint64_t>> sketches;
    for (uint32_t u = 0; u < n_phys; ++u) {
        if (result.phasing_nodes.count(u)) continue;
        const std::string& seq = graph.get_sequence(u);
        if (static_cast<int>(seq.size()) < k) continue;
        auto sk = compute_sketch(seq, k, sketch_size);
        if (!sk.empty()) sketches.emplace(u, std::move(sk));
    }

    std::unordered_map<uint32_t, std::pair<uint32_t, double>> best;
    for (const auto& kv : sketches) {
        const uint32_t u = kv.first;
        const auto& sk_u = kv.second;

        std::unordered_set<uint32_t> siblings;
        for (uint32_t p : in_phys[u])
            for (uint32_t c : out_phys[p])
                if (c != u) siblings.insert(c);
        for (uint32_t s : out_phys[u])
            for (uint32_t pp : in_phys[s])
                if (pp != u) siblings.insert(pp);

        uint32_t best_v = u;
        double   best_j = -1.0;
        for (uint32_t v : siblings) {
            auto it = sketches.find(v);
            if (it == sketches.end()) continue;
            double j = estimate_jaccard(sk_u, it->second);
            if (j >= jaccard_threshold && j > best_j) { best_j = j; best_v = v; }
        }
        if (best_v != u) best[u] = {best_v, best_j};
    }

    for (const auto& kv : best) {
        const uint32_t u = kv.first;
        const uint32_t v = kv.second.first;
        if (u >= v) continue;
        auto it = best.find(v);
        if (it != best.end() && it->second.first == u) {
            result.phasing_nodes.insert(u);
            result.phasing_nodes.insert(v);
            result.alt_map[u] = v;
            result.alt_map[v] = u;
        }
    }

    return result;
}

BubbleDetector::AssemblerStyle
BubbleDetector::detect_style(const Graph& graph) const {
    const std::regex re_classic  {R"(^\d+[FR](_\d+)?$)"};       // Shasta classic
    const std::regex re_pr       {R"(^PR\.\d+\.\d+$)"};         // Shasta PR-style
    const std::regex re_dot      {R"(^.+\.\d+$)"};              // generic dot-suffix
    const std::regex re_hifiasm  {R"(^h[12]tg\d+l$)"};          // Hifiasm: h1tg000001l / h2tg000001l
    const std::regex re_verkko   {R"(^haplotype[12]-\d+$)"};    // Verkko:  haplotype1-0000001 / haplotype2-0000001

    int votes[5] = {0, 0, 0, 0, 0};  // Classic, PRStyle, DotSuffix, HifiAsm, Verkko
    size_t limit = std::min(graph.get_num_nodes(), size_t{1000});

    for (size_t i = 0; i < limit; ++i) {
        const std::string& name = graph.get_name(static_cast<uint32_t>(i));
        if (std::regex_match(name, re_classic)) ++votes[0];
        if (std::regex_match(name, re_pr))      ++votes[1];
        if (std::regex_match(name, re_dot))     ++votes[2];
        if (std::regex_match(name, re_hifiasm)) ++votes[3];
        if (std::regex_match(name, re_verkko))  ++votes[4];
    }

    int best = static_cast<int>(
        std::max_element(votes, votes + 5) - votes);

    switch (best) {
        case 0:  return AssemblerStyle::Classic;
        case 1:  return AssemblerStyle::PRStyle;
        case 2:  return AssemblerStyle::DotSuffix;
        case 3:  return AssemblerStyle::HifiAsmStyle;
        case 4:  return AssemblerStyle::VerkkoStyle;
        default: return AssemblerStyle::DotSuffix;
    }
}

std::string BubbleDetector::base_name(const std::string& name,
                                       AssemblerStyle style) const {
    switch (style) {
        case AssemblerStyle::Classic: {
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
        case AssemblerStyle::PRStyle:
        case AssemblerStyle::DotSuffix: {
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
        case AssemblerStyle::HifiAsmStyle:
            // TODO: Strip h1/h2 prefix (e.g. "h1tg000001l" → "tg000001l").

            return {};
        case AssemblerStyle::VerkkoStyle:
            // TODO: Strip haplotype[12]- prefix (e.g. "haplotype1-0000001" → "0000001").
 
            return {};
    }
    return {};
}

BubbleResult BubbleDetector::get_alts_from_shasta_names(const Graph& graph) const {
    AssemblerStyle style = detect_style(graph);

    std::unordered_map<std::string, std::vector<uint32_t>> groups;
    for (uint32_t id = 0; id < graph.get_num_nodes(); ++id) {
        std::string base = base_name(graph.get_name(id), style);
        if (!base.empty())
            groups[base].push_back(id);
    }

    BubbleResult result;
    for (auto& [base, ids] : groups) {
        if (ids.size() != 2) continue;  

        uint32_t a = ids[0], b = ids[1];
        result.phasing_nodes.insert(a);
        result.phasing_nodes.insert(b);
        result.alt_map[a] = b;
        result.alt_map[b] = a;
    }

    return result;
}

namespace {

// Find the index of a header column by exact name. Returns -1 if not found.
int find_col(const std::vector<std::string>& header, const std::string& name) {
    for (size_t i = 0; i < header.size(); ++i)
        if (header[i] == name) return static_cast<int>(i);
    return -1;
}

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == ',') { out.push_back(std::move(cur)); cur.clear(); }
        else cur.push_back(c);
    }
    out.push_back(std::move(cur));
    return out;
}

} // namespace

BubbleResult BubbleDetector::from_shasta_csv(const Graph& graph,
                                              const std::string& csv_path) const {
    std::ifstream in(csv_path);
    if (!in) throw std::runtime_error("Cannot open Shasta CSV: " + csv_path);

    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("Empty CSV: " + csv_path);

    auto header = split_csv(line);
    const int col_name      = find_col(header, "Name");
    const int col_ploidy    = find_col(header, "Ploidy");
    const int col_chain     = find_col(header, "Bubble chain");
    const int col_pos       = find_col(header, "Position in bubble chain");
    const int col_component = find_col(header, "Component");
    const int col_haplotype = find_col(header, "Haplotype");

    if (col_name < 0 || col_ploidy < 0 || col_chain < 0 ||
        col_pos  < 0 || col_component < 0 || col_haplotype < 0) {
        throw std::runtime_error(
            "Shasta CSV missing required columns "
            "(Name, Ploidy, Bubble chain, Position in bubble chain, Component, Haplotype): "
            + csv_path);
    }

    // (chain, position, component) -> {hap0_id, hap1_id} (UINT32_MAX = unset)
    using Key = std::tuple<std::string, std::string, std::string>;
    struct KeyHash {
        size_t operator()(const Key& k) const {
            size_t h = std::hash<std::string>{}(std::get<0>(k));
            h ^= std::hash<std::string>{}(std::get<1>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<std::string>{}(std::get<2>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };
    std::unordered_map<Key, std::pair<uint32_t,uint32_t>, KeyHash> groups;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto cols = split_csv(line);
        if (static_cast<int>(cols.size()) <= col_haplotype) continue;
        if (cols[col_ploidy] != "2") continue;

        const std::string& name = cols[col_name];
        auto it = std::find_if(name.begin(), name.end(),
                               [](unsigned char c) { return !std::isspace(c); });
        if (it == name.end()) continue;

        uint32_t id;
        try { id = graph.get_id(name); }
        catch (...) { continue; }  // segment not in graph (e.g. filtered out)

        Key k(cols[col_chain], cols[col_pos], cols[col_component]);
        auto ins = groups.emplace(k, std::make_pair(UINT32_MAX, UINT32_MAX));
        auto& slot = ins.first->second;
        if      (cols[col_haplotype] == "0") slot.first  = id;
        else if (cols[col_haplotype] == "1") slot.second = id;
    }

    BubbleResult result;
    for (const auto& kv : groups) {
        uint32_t a = kv.second.first;
        uint32_t b = kv.second.second;
        if (a == UINT32_MAX || b == UINT32_MAX) continue;  // incomplete pair
        result.phasing_nodes.insert(a);
        result.phasing_nodes.insert(b);
        result.alt_map[a] = b;
        result.alt_map[b] = a;
    }
    return result;
}
