#include "Graph.hpp"
#include "ContactMatrix.hpp"
#include "BubbleDetector.hpp"
#include "TrioBinner.hpp"
#include "Phaser.hpp"
#include "Chainer.hpp"
#include "HamiltonianChainer.hpp"
#include "UnzippedGraph.hpp"
#include "ResultWriter.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread_pool/thread_pool.hpp>

namespace {
void print_usage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name
        << " -g <file.gfa> -b <file.bam> [options]\n"
        << "\n"
        << "Required:\n"
        << "  -g, --gfa              <path>           Assembly graph (GFA format)\n"
        << "  -b, --bam              <path>           Hi-C / Pore-C alignments (BAM format)\n"
        << "\n"
        << "Optional:\n"
        << "  -f, --fasta            <path>           Companion FASTA (needed when GFA uses '*')\n"
        << "  -p, --pat-yak          <path>           Paternal k-mer database (.yak) for trio binning\n"
        << "  -m, --mat-yak          <path>           Maternal k-mer database (.yak) for trio binning\n"
        << "      --bubble-mode      shasta|minhash|csv  Bubble detection mode (default: shasta)\n"
        << "      --shasta-csv       <path>           Shasta Assembly-Phased.csv (required if --bubble-mode csv)\n"
        << "      --kmer-size        <int>            k-mer size for MinHash (default: 16)\n"
        << "      --sketch-size      <int>            MinHash sketch size (default: 1000)\n"
        << "      --jaccard-threshold <float>         MinHash Jaccard threshold (default: 0.2)\n"
        << "  -t, --threads          <int>            Number of threads (default: 4)\n"
        << "      --out-dir          <path>           Output directory (default: \".\")\n"
        << "      --prefix           <name>           Output filename prefix (default: GFA basename)\n"
        << "      --simple-chainer                    Use greedy Chainer (default: HamiltonianChainer)\n"
        << "      --hamiltonian                       (kept for backward compat; HamiltonianChainer is default)\n"
        << "  -h, --help                              Show this message\n";
}
} // namespace

int main(int argc, char* argv[]) {
    std::string gfa_path;
    std::string bam_path;
    std::string fasta_path;
    std::string pat_yak_path;
    std::string mat_yak_path;
    std::string bubble_mode       = "shasta";
    std::string shasta_csv_path;
    int         kmer_size         = 16;
    int         sketch_size       = 1000;
    double      jaccard_threshold = 0.2;
    int         num_threads       = 4;
    std::string out_dir           = ".";
    std::string prefix;
    bool        use_simple_chainer = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }

        auto require_next = [&](const std::string& flag) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << flag << "\n";
                print_usage(argv[0]);
                std::exit(1);
            }
            return argv[++i];
        };

        if (arg == "-g" || arg == "--gfa")     { gfa_path     = require_next(arg); continue; }
        if (arg == "-b" || arg == "--bam")     { bam_path     = require_next(arg); continue; }
        if (arg == "-f" || arg == "--fasta")   { fasta_path   = require_next(arg); continue; }
        if (arg == "-p" || arg == "--pat-yak") { pat_yak_path = require_next(arg); continue; }
        if (arg == "-m" || arg == "--mat-yak") { mat_yak_path = require_next(arg); continue; }

        if (arg == "--bubble-mode") {
            bubble_mode = require_next(arg);
            if (bubble_mode != "shasta" && bubble_mode != "minhash" && bubble_mode != "csv") {
                std::cerr << "--bubble-mode must be 'shasta', 'minhash', or 'csv'\n";
                return 1;
            }
            continue;
        }
        if (arg == "--shasta-csv") { shasta_csv_path = require_next(arg); continue; }
        if (arg == "-t" || arg == "--threads") { num_threads = std::stoi(require_next(arg)); continue; }
        if (arg == "--kmer-size")          { kmer_size         = std::stoi(require_next(arg)); continue; }
        if (arg == "--sketch-size")        { sketch_size       = std::stoi(require_next(arg)); continue; }
        if (arg == "--jaccard-threshold")  { jaccard_threshold = std::stod(require_next(arg)); continue; }
        if (arg == "--out-dir")            { out_dir           = require_next(arg); continue; }
        if (arg == "--prefix")             { prefix            = require_next(arg); continue; }
        if (arg == "--hamiltonian")        { /* default; kept for backward compat */ continue; }
        if (arg == "--simple-chainer")     { use_simple_chainer = true; continue; }

        std::cerr << "Unknown argument: " << arg << "\n";
        print_usage(argv[0]);
        return 1;
    }

    if (gfa_path.empty() || bam_path.empty()) {
        std::cerr << "Both --gfa and --bam are required.\n";
        print_usage(argv[0]);
        return 1;
    }
    if (pat_yak_path.empty() != mat_yak_path.empty()) {
        std::cerr << "Both --pat-yak and --mat-yak must be provided together.\n";
        return 1;
    }
    if (bubble_mode == "csv" && shasta_csv_path.empty()) {
        std::cerr << "--bubble-mode csv requires --shasta-csv <path>.\n";
        return 1;
    }

    auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

    try {
        // --- Stage 1: parse graph ---
        Graph graph;
        graph.load_from_gfa(gfa_path);
        std::cout << "Loaded graph: " << graph.get_num_nodes() << " nodes\n";

        if (!fasta_path.empty()) {
            graph.load_fasta(fasta_path);
            std::cout << "Loaded sequences from FASTA: " << fasta_path << "\n";
        }

        // --- Stage 1: build contact matrix ---
        ContactMatrix contacts;
        contacts.load_contacts(bam_path, graph, num_threads);
        contacts.build_csr(graph.get_num_nodes(), thread_pool);
        std::cout << "Built contact matrix.\n";

        // --- Stage 1b (optional): trio binning — computed now, applied after bubble detection ---
        std::vector<TrioScores> trio_scores;
        if (!pat_yak_path.empty()) {
            TrioBinner binner(pat_yak_path, mat_yak_path);
            trio_scores = binner.compute_scores(graph, thread_pool);
            uint64_t total_pat = 0, total_mat = 0;
            for (const auto& s : trio_scores) { total_pat += s.pat_count; total_mat += s.mat_count; }
            std::cout << "Trio binning: paternal k-mers = " << total_pat
                      << ", maternal = " << total_mat << "\n";
        }

        // --- Stage 2: bubble / homology detection ---
        if (bubble_mode == "minhash") {
            bool has_seq = false;
            for (size_t i = 0; i < graph.get_num_nodes(); ++i)
                if (!graph.get_sequence(static_cast<uint32_t>(i)).empty()) { has_seq = true; break; }
            if (!has_seq)
                throw std::runtime_error(
                    "--bubble-mode minhash requires sequences. "
                    "Provide a companion FASTA with -f/--fasta.");
        }

        BubbleDetector detector;
        BubbleResult bubbles;
        if      (bubble_mode == "shasta")  bubbles = detector.get_alts_from_shasta_names(graph);
        else if (bubble_mode == "csv")     bubbles = detector.from_shasta_csv(graph, shasta_csv_path);
        else                                bubbles = detector.find_unlabeled_alts(
                                                graph, kmer_size, sketch_size, jaccard_threshold, thread_pool);

        std::cout << "Bubble detection (" << bubble_mode << "): "
                  << bubbles.phasing_nodes.size() << " phasing nodes, "
                  << bubbles.alt_map.size() / 2   << " bubble pairs\n";

        contacts.filter_to_phasing_nodes(bubbles.phasing_nodes, graph.get_num_nodes());
        std::cout << "Contact matrix filtered to bubble nodes.\n";

        // --- Stage 2b (optional): lock trio-dominant phasing nodes ---
        if (!trio_scores.empty()) {
            size_t locked = apply_trio_constraints(graph, bubbles.phasing_nodes, trio_scores);
            std::cout << "Trio constraints locked " << locked << " / "
                      << bubbles.phasing_nodes.size() << " phasing nodes.\n";
        }

        // --- Stage 3: Monte Carlo phasing ---
        monte_carlo_phase(graph, contacts, bubbles, PhaserConfig{}, thread_pool);
        std::cout << "Monte Carlo phasing complete.\n";

        // --- Stage 4: chaining + unzipping ---
        ChainResult chains = use_simple_chainer
            ? Chainer().generate_chain_paths(graph, bubbles, &contacts)
            : HamiltonianChainer(thread_pool).generate_chain_paths(graph, bubbles, contacts);
        std::cout << "Chained: ph0=" << chains.phase_0.size()
                  << " ph1=" << chains.phase_1.size()
                  << " unphased=" << chains.unphased.size() << "\n";

        UnzippedGraph uz = UnzippedGraph::build(graph, bubbles, chains);
        std::cout << "Unzipped graph: " << uz.segments().size() << " segments, "
                  << uz.links().size() << " links\n";

        // --- Stage 5: serialize results ---
        if (prefix.empty()) {
            prefix = std::filesystem::path(gfa_path).stem().string();
            if (prefix.empty()) prefix = "out";
        }
        OutputConfig out_cfg{out_dir, prefix};
        ResultWriter writer(thread_pool);
        writer.write_chained_gfa (graph, chains, out_cfg);
        writer.write_unzipped_gfa(graph, uz,     out_cfg);
        writer.write_fastas      (graph, uz,     out_cfg);
        std::cout << "Wrote results under " << out_dir << " with prefix '" << prefix << "'\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
