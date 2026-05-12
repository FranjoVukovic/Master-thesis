#include "ResultWriter.hpp"

#include "Graph.hpp"

#include <filesystem>
#include <fstream>
#include <future>
#include <stdexcept>
#include <unordered_set>

namespace {

std::string join_path(const std::string& dir, const std::string& name) {
    if (dir.empty() || dir == ".") return name;
    return (std::filesystem::path(dir) / name).string();
}

void ensure_dir(const std::string& dir) {
    if (dir.empty() || dir == ".") return;
    std::filesystem::create_directories(dir);
}

const char* strand(bool rev) { return rev ? "-" : "+"; }

void write_fasta_record(std::ofstream& out,
                        const std::string& name,
                        const std::string& seq) {
    out << ">" << name << "\n";
    const size_t line = 60;
    for (size_t i = 0; i < seq.size(); i += line) {
        out.write(seq.data() + i, std::min(line, seq.size() - i));
        out << "\n";
    }
}

} // namespace

ResultWriter::ResultWriter(std::shared_ptr<thread_pool::ThreadPool> pool)
    : pool_(std::move(pool)) {}

void ResultWriter::write_chained_gfa(const Graph& graph,
                                     const ChainResult& chains,
                                     const OutputConfig& cfg) const {
    ensure_dir(cfg.out_dir);
    const std::string path = join_path(cfg.out_dir, cfg.prefix + ".chained.gfa");
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open output: " + path);

    const size_t N = graph.get_num_nodes();

    // S-lines.
    for (uint32_t i = 0; i < N; ++i) {
        const std::string& seq = graph.get_sequence(i);
        out << "S\t" << graph.get_name(i) << "\t"
            << (seq.empty() ? "*" : seq);
        const uint32_t len = graph.get_length(i)
                             ? graph.get_length(i)
                             : static_cast<uint32_t>(seq.size());
        if (len > 0) out << "\tLN:i:" << len;
        if (graph.get_coverage(i) > 0) out << "\tRC:i:" << graph.get_coverage(i);
        out << "\tPH:i:" << static_cast<int>(graph.get_phase(i));
        out << "\n";
    }

    // L-lines: emit once per unordered pair using the edge stored at min(u,v).
    std::unordered_set<uint64_t> emitted;
    auto key = [](uint32_t a, uint32_t b) {
        return (static_cast<uint64_t>(a) << 32) | b;
    };
    for (uint32_t u = 0; u < N; ++u) {
        auto range = graph.get_neighbors(u);
        for (auto it = range.first; it != range.second; ++it) {
            uint32_t v = it->target_id;
            uint32_t a = std::min(u, v), b = std::max(u, v);
            if (!emitted.insert(key(a, b)).second) continue;
            if (u != a) continue;  // emit only when we are at the lower id
            out << "L\t" << graph.get_name(u) << "\t" << strand(it->source_rev)
                << "\t"  << graph.get_name(v) << "\t" << strand(it->target_rev)
                << "\t"  << it->overlap_length << "M\n";
        }
    }

    // P-lines: one walk per chain.
    auto emit_paths = [&](const std::vector<Chain>& cs, const std::string& tag) {
        for (size_t i = 0; i < cs.size(); ++i) {
            const Chain& c = cs[i];
            if (c.nodes.empty()) continue;
            out << "P\t" << tag << "_" << i << "\t";
            for (size_t j = 0; j < c.nodes.size(); ++j) {
                if (j) out << ",";
                out << graph.get_name(c.nodes[j]) << (c.reverse[j] ? "-" : "+");
            }
            out << "\t*\n";
        }
    };
    emit_paths(chains.phase_0,  "chain_ph0");
    emit_paths(chains.phase_1,  "chain_ph1");
    emit_paths(chains.unphased, "chain_un");
}

void ResultWriter::write_unzipped_gfa(const Graph& graph,
                                      const UnzippedGraph& uz,
                                      const OutputConfig& cfg) const {
    ensure_dir(cfg.out_dir);
    const std::string path = join_path(cfg.out_dir, cfg.prefix + ".unzipped.gfa");
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot open output: " + path);

    for (const auto& s : uz.segments()) {
        const std::string& seq = graph.get_sequence(s.src_id);
        out << "S\t" << s.name << "\t" << (seq.empty() ? "*" : seq);
        if (s.length > 0)   out << "\tLN:i:" << s.length;
        if (s.coverage > 0) out << "\tRC:i:" << s.coverage;
        out << "\tPH:i:" << static_cast<int>(s.phase) << "\n";
    }
    for (const auto& l : uz.links()) {
        out << "L\t" << uz.segments()[l.source_idx].name << "\t" << strand(l.source_rev)
            << "\t"  << uz.segments()[l.target_idx].name << "\t" << strand(l.target_rev)
            << "\t"  << l.overlap_length << "M\n";
    }
}

void ResultWriter::write_fastas(const Graph& graph,
                                const UnzippedGraph& uz,
                                const OutputConfig& cfg) const {
    ensure_dir(cfg.out_dir);
    const std::string p0 = join_path(cfg.out_dir, cfg.prefix + ".phase_0.fasta");
    const std::string p1 = join_path(cfg.out_dir, cfg.prefix + ".phase_1.fasta");
    const std::string pu = join_path(cfg.out_dir, cfg.prefix + ".unphased.fasta");

    auto write_partition = [&](const std::string& path, int8_t phase) {
        std::ofstream out(path);
        if (!out) throw std::runtime_error("Cannot open output: " + path);
        for (const auto& s : uz.segments()) {
            if (s.phase != phase) continue;
            write_fasta_record(out, s.name, graph.get_sequence(s.src_id));
        }
    };

    if (pool_ && pool_->num_threads() >= 2) {
        auto f0 = pool_->Submit([&]() { write_partition(p0,  0); return 0; });
        auto f1 = pool_->Submit([&]() { write_partition(p1,  1); return 0; });
        auto fu = pool_->Submit([&]() { write_partition(pu, -1); return 0; });
        f0.get(); f1.get(); fu.get();
    } else {
        write_partition(p0,  0);
        write_partition(p1,  1);
        write_partition(pu, -1);
    }
}
