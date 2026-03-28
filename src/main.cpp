#include "Graph.hpp"
#include "ContactMatrix.hpp"
#include "TrioBinner.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void print_usage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name
        << " -g <file.gfa> -b <file.bam> [-p <pat.yak> -m <mat.yak>]\n"
        << "\n"
        << "Required:\n"
        << "  -g, --gfa   <path>   Assembly graph (GFA format)\n"
        << "  -b, --bam   <path>   Hi-C / Pore-C alignments (BAM format)\n"
        << "\n"
        << "Optional:\n"
        << "  -p, --pat-yak <path> Paternal k-mer database (yak format) for trio binning\n"
        << "  -m, --mat-yak <path> Maternal k-mer database (yak format) for trio binning\n"
        << "  -h, --help           Show this message\n";
}
} // namespace

int main(int argc, char* argv[]) {
    std::string gfa_path;
    std::string bam_path;
    std::string pat_yak_path;
    std::string mat_yak_path;

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

        if (arg == "-g" || arg == "--gfa") { gfa_path     = require_next(arg); continue; }
        if (arg == "-b" || arg == "--bam") { bam_path     = require_next(arg); continue; }
        if (arg == "-p" || arg == "--pat-yak") { pat_yak_path = require_next(arg); continue; }
        if (arg == "-m" || arg == "--mat-yak") { mat_yak_path = require_next(arg); continue; }

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

    try {
        // Stage 1: parse graph and contact matrix.
        Graph graph;
        graph.load_from_gfa(gfa_path);
        std::cout << "Loaded graph: " << graph.get_num_nodes() << " nodes\n";

        ContactMatrix contacts;
        contacts.load_contacts(bam_path, graph);
        contacts.build_csr(graph.get_num_nodes());
        std::cout << "Built contact matrix.\n";

        // Stage 1b (optional): trio binning scores.
        if (!pat_yak_path.empty()) {
            TrioBinner binner(pat_yak_path, mat_yak_path);
            auto scores = binner.compute_scores(graph);

            uint64_t total_pat = 0, total_mat = 0;
            for (const auto& s : scores) {
                total_pat += s.pat_count;
                total_mat += s.mat_count;
            }
            std::cout << "Trio binning: total paternal k-mers = " << total_pat
                      << ", maternal = " << total_mat << "\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
