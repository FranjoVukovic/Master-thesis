#include <gtest/gtest.h>

#include "ContactMatrix.hpp"
#include "Graph.hpp"

#include <cstdint>
#include <string>

TEST(GraphTest, GfaParsingAndStructure) {
    Graph graph;
    graph.load_from_gfa(std::string(TEST_DATA_DIR) + "/test.gfa");

    const size_t num_nodes = graph.get_num_nodes();
    ASSERT_EQ(num_nodes, 3u);

    EXPECT_EQ(graph.get_name(0), "contig1");
    EXPECT_EQ(graph.get_name(1), "contig2");
    EXPECT_EQ(graph.get_name(2), "contig3");

    const auto neighbors0 = graph.get_neighbors(0);
    bool edge_found = false;
    for (auto it = neighbors0.first; it != neighbors0.second; ++it) {
        if (it->target_id == 1 && it->overlap_length == 0) {
            edge_found = true;
            break;
        }
    }
    EXPECT_TRUE(edge_found);
}

TEST(ContactMatrixTest, BamParsingAndCsrConstruction) {
    Graph graph;
    graph.load_from_gfa(std::string(TEST_DATA_DIR) + "/test.gfa");

    ContactMatrix contacts;
    contacts.load_contacts(std::string(TEST_DATA_DIR) + "/test.bam", graph);
    
    const size_t num_nodes = graph.get_num_nodes();
    contacts.build_csr(num_nodes);

    EXPECT_EQ(contacts.get_contact(0, 1), 1u);
    EXPECT_EQ(contacts.get_contact(1, 2), 1u);
    EXPECT_EQ(contacts.get_contact(0, 2), 0u);
}