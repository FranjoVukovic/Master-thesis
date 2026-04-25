#include "TrioBinner.hpp"
#include "Graph.hpp"

// yak C API (compiled from FetchContent source)
extern "C" {
#include "yak.h"
}

#include <future>
#include <stdexcept>
#include <string>

// ---------------------------------------------------------------------------
// Canonical k-mer hashing — mirrors yak's internal logic.
// yak stores the minimum of the forward and reverse-complement k-mer hashes.
// We use the same hash mix as yak (hash64 from yak.h / sketch.c).
// ---------------------------------------------------------------------------

static inline uint64_t yak_hash64_local(uint64_t key, uint64_t mask) {
    key = (~key + (key << 21)) & mask;
    key ^= key >> 24;
    key  = ((key + (key << 3)) + (key << 8)) & mask;
    key ^= key >> 14;
    key  = ((key + (key << 2)) + (key << 4)) & mask;
    key ^= key >> 28;
    key  = (key + (key << 31)) & mask;
    return key;
}

// Nucleotide → 2-bit encoding used by yak (A=0, C=1, G=2, T=3).
static const uint8_t nt_table[256] = {
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 0-15
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 16-31
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 32-47
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 48-63
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,  // 64-79  (A,C,G)
    4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4,  // 80-95  (T)
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,  // 96-111 (a,c,g)
    4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4,  // 112-127 (t)
};

// ---------------------------------------------------------------------------

TrioBinner::TrioBinner(const std::string& pat_yak, const std::string& mat_yak, int k)
    : pat_db_(nullptr), mat_db_(nullptr), k_(k)
{
    pat_db_ = static_cast<void*>(yak_ch_restore(pat_yak.c_str()));
    if (!pat_db_)
        throw std::runtime_error("Failed to open paternal yak database: " + pat_yak);

    mat_db_ = static_cast<void*>(yak_ch_restore(mat_yak.c_str()));
    if (!mat_db_) {
        yak_ch_destroy(static_cast<yak_ch_t*>(pat_db_));
        throw std::runtime_error("Failed to open maternal yak database: " + mat_yak);
    }
}

TrioBinner::~TrioBinner() {
    if (pat_db_) yak_ch_destroy(static_cast<yak_ch_t*>(pat_db_));
    if (mat_db_) yak_ch_destroy(static_cast<yak_ch_t*>(mat_db_));
}

std::vector<TrioScores> TrioBinner::compute_scores(const Graph& graph,
        std::shared_ptr<thread_pool::ThreadPool> pool) const {
    const size_t n = graph.get_num_nodes();
    std::vector<TrioScores> scores(n);

    if (!pool) {
        for (size_t i = 0; i < n; ++i) {
            const std::string& seq = graph.get_sequence(static_cast<uint32_t>(i));
            if (seq.empty()) continue;
            scores[i].pat_count = count_matching_kmers(seq, pat_db_);
            scores[i].mat_count = count_matching_kmers(seq, mat_db_);
        }
        return scores;
    }

    const size_t num_threads = pool->num_threads();
    const size_t chunk_size = std::max<size_t>(1, (n + num_threads - 1) / num_threads);
    std::vector<std::future<void>> futures;

    for (size_t start = 0; start < n; start += chunk_size) {
        size_t end = std::min(start + chunk_size, n);
        futures.emplace_back(pool->Submit(
            [this, &graph, &scores](size_t begin, size_t end) {
                for (size_t i = begin; i < end; ++i) {
                    const std::string& seq = graph.get_sequence(static_cast<uint32_t>(i));
                    if (seq.empty()) continue;
                    scores[i].pat_count = count_matching_kmers(seq, pat_db_);
                    scores[i].mat_count = count_matching_kmers(seq, mat_db_);
                }
            },
            start, end));
    }
    for (auto& f : futures) f.get();

    return scores;
}

uint32_t TrioBinner::count_matching_kmers(const std::string& seq, void* db) const {
    yak_ch_t* ch = static_cast<yak_ch_t*>(db);
    const uint64_t mask = (1ULL << (2 * k_)) - 1;

    uint64_t fwd = 0, rev = 0;
    int valid_bases = 0;
    uint32_t count = 0;

    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t b = nt_table[static_cast<unsigned char>(seq[i])];
        if (b >= 4) {
            // Ambiguous base — reset the k-mer window.
            fwd = rev = 0;
            valid_bases = 0;
            continue;
        }

        // Shift new base into forward and reverse-complement accumulators.
        fwd = ((fwd << 2) | b) & mask;
        rev = ((rev >> 2) | (static_cast<uint64_t>(3 ^ b) << (2 * (k_ - 1)))) & mask;
        ++valid_bases;

        if (valid_bases < k_) continue;

        // Canonical k-mer is the hash of whichever strand gives the smaller hash.
        uint64_t hf = yak_hash64_local(fwd, mask);
        uint64_t hr = yak_hash64_local(rev, mask);
        uint64_t canonical = (hf < hr) ? hf : hr;

        // Query the yak hash table. yak_ch_get returns the count (0 if absent).
        if (yak_ch_get(ch, canonical) > 0)
            ++count;
    }

    return count;
}
