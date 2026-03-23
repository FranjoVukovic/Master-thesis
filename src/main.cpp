#include "Graph.hpp"
#include "ContactMatrix.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void print_usage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name << " --gfa <path/to/file.gfa> --bam <path/to/file.bam>\n"
        << "   or: " << program_name << " -g <path/to/file.gfa> -b <path/to/file.bam>\n";
}
}

int main(int argc, char* argv[]) {
    std::string gfa_path;
    std::string bam_path;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }

        if (arg == "-g" || arg == "--gfa") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << arg << "\n";
                print_usage(argv[0]);
                return 1;
            }
            gfa_path = argv[++i];
            continue;
        }

        if (arg == "-b" || arg == "--bam") {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << arg << "\n";
                print_usage(argv[0]);
                return 1;
            }
            bam_path = argv[++i];
            continue;
        }

        std::cerr << "Unknown argument: " << arg << "\n";
        print_usage(argv[0]);
        return 1;
    }

    if (gfa_path.empty() || bam_path.empty()) {
        std::cerr << "Both --gfa and --bam are required.\n";
        print_usage(argv[0]);
        return 1;
    }

    try {
        Graph graph;
        graph.load_from_gfa(gfa_path);

        ContactMatrix contacts;
        contacts.load_contacts(bam_path, graph);

        size_t num_nodes = graph.get_num_nodes();
        contacts.build_csr(num_nodes);

        std::cout << "Successfully loaded GFA: " << gfa_path << "\n";
        std::cout << "Successfully loaded BAM: " << bam_path << "\n";

        if (num_nodes == 3) {
            std::cout << "\n--- Structure Verification (Dummy Data Test) ---\n";

            // 1. Verify Node Names
            std::string name0 = graph.get_name(0); // Expected: contig1
            std::string name1 = graph.get_name(1); // Expected: contig2
            std::string name2 = graph.get_name(2); // Expected: contig3
            
            std::cout << "Nodes loaded: " << name0 << ", " << name1 << ", " << name2 << "\n";

            // 2. Verify Graph Edges (contig1 -> contig2)
            auto neighbors0 = graph.get_neighbors(0);
            bool edge_found = false;
            for (auto it = neighbors0.first; it != neighbors0.second; ++it) {
                if (it->target_id == 1 && it->overlap_length == 0) {
                    edge_found = true;
                    break;
                }
            }
            std::cout << "Edge " << name0 << " -> " << name1 << " present: " 
                      << (edge_found ? "Yes" : "No") << "\n";

            // 3. Verify Contact Matrix values
            uint32_t contacts_0_1 = contacts.get_contact(0, 1);
            uint32_t contacts_1_2 = contacts.get_contact(1, 2);
            uint32_t contacts_0_2 = contacts.get_contact(0, 2); // Should be 0

            std::cout << "Contacts " << name0 << "-" << name1 << " (Expected 1): " << contacts_0_1 << "\n";
            std::cout << "Contacts " << name1 << "-" << name2 << " (Expected 1): " << contacts_1_2 << "\n";
            std::cout << "Contacts " << name0 << "-" << name2 << " (Expected 0): " << contacts_0_2 << "\n";
            
            // 4. Final Assessment
            if (contacts_0_1 == 1 && contacts_1_2 == 1 && contacts_0_2 == 0 && edge_found) {
                std::cout << "-> ALL TESTS PASSED.\n";
            } else {
                std::cout << "-> SOME TESTS FAILED.\n";
            }
        } else if (num_nodes > 0) {
            std::cout << "\n[Note: Loaded " << num_nodes 
                      << " nodes. Dummy data test skipped as it requires exactly 3 nodes.]\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}