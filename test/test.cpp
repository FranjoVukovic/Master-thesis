#include <gtest/gtest.h>

#include "BubbleDetector.hpp"
#include "ContactMatrix.hpp"
#include "Graph.hpp"
#include "Phaser.hpp"
#include "TrioBinner.hpp"

#include <cstdint>
#include <string>
#include <unordered_set>

class UnifiedTest : public ::testing::Test {
protected:
    Graph          graph;
    BubbleDetector det;

    void SetUp() override {
        graph.load_from_gfa(std::string(TEST_DATA_DIR) + "/test_unified.gfa");
    }
};



TEST_F(UnifiedTest, GfaParsingAndLinearChain) {
    EXPECT_GT(graph.get_num_nodes(), 0u);

    const uint32_t lc_a = graph.get_id("lc_a");
    const uint32_t lc_b = graph.get_id("lc_b");
    const uint32_t lc_c = graph.get_id("lc_c");
    const uint32_t lc_d = graph.get_id("lc_d");

    EXPECT_EQ(graph.get_name(lc_a), "lc_a");
    EXPECT_EQ(graph.get_name(lc_b), "lc_b");
    EXPECT_EQ(graph.get_name(lc_c), "lc_c");
    EXPECT_EQ(graph.get_name(lc_d), "lc_d");

    auto has_neighbor = [&](uint32_t u, uint32_t v) {
        auto range = graph.get_neighbors(u);
        for (auto it = range.first; it != range.second; ++it) {
            if (it->target_id == v && it->overlap_length == 0) return true;
        }
        return false;
    };

    EXPECT_TRUE(has_neighbor(lc_a, lc_b));
    EXPECT_TRUE(has_neighbor(lc_b, lc_c));
    EXPECT_TRUE(has_neighbor(lc_c, lc_d));

    EXPECT_TRUE(has_neighbor(lc_b, lc_a));
    EXPECT_TRUE(has_neighbor(lc_c, lc_b));
    EXPECT_TRUE(has_neighbor(lc_d, lc_c));
}


TEST_F(UnifiedTest, BamParsingAndCsrConstruction) {
    ContactMatrix contacts;
    contacts.load_contacts(std::string(TEST_DATA_DIR) + "/test_unified.bam", graph);
    contacts.build_csr(graph.get_num_nodes());

    const uint32_t lc_a = graph.get_id("lc_a");
    const uint32_t lc_b = graph.get_id("lc_b");
    const uint32_t lc_c = graph.get_id("lc_c");

    EXPECT_EQ(contacts.get_contact(lc_a, lc_b), 1u);
    EXPECT_EQ(contacts.get_contact(lc_b, lc_c), 1u);
    EXPECT_EQ(contacts.get_contact(lc_a, lc_c), 0u);
}

class UnifiedBubbleTest : public UnifiedTest {
protected:
    BubbleResult r;

    void SetUp() override {
        UnifiedTest::SetUp();
        r = det.find_unlabeled_alts(graph, 16, 1000, 0.0);
    }
};

TEST_F(UnifiedBubbleTest, SNPBubble) {
    const uint32_t snp_A = graph.get_id("snp_A");
    const uint32_t snp_B = graph.get_id("snp_B");

    EXPECT_EQ(r.phasing_nodes.count(snp_A), 1u);
    EXPECT_EQ(r.phasing_nodes.count(snp_B), 1u);
    EXPECT_EQ(r.alt_map.at(snp_A), snp_B);
    EXPECT_EQ(r.alt_map.at(snp_B), snp_A);
}

TEST_F(UnifiedBubbleTest, MultiNodeSuperbubble) {
    const uint32_t mb_A = graph.get_id("mb_A");
    const uint32_t mb_B = graph.get_id("mb_B");
    const uint32_t mb_C = graph.get_id("mb_C");
    const uint32_t mb_D = graph.get_id("mb_D");

    EXPECT_EQ(r.phasing_nodes.count(mb_A), 1u);
    EXPECT_EQ(r.phasing_nodes.count(mb_B), 1u);
    EXPECT_EQ(r.phasing_nodes.count(mb_C), 1u);
    EXPECT_EQ(r.phasing_nodes.count(mb_D), 1u);

    EXPECT_EQ(r.alt_map.at(mb_A), mb_B);
    EXPECT_EQ(r.alt_map.at(mb_C), mb_B);
    EXPECT_EQ(r.alt_map.at(mb_B), mb_A);
    EXPECT_EQ(r.alt_map.at(mb_D), mb_A);
}

TEST_F(UnifiedBubbleTest, ThreeBranchComplexSuperbubble) {
    const uint32_t tb_A = graph.get_id("tb_A");
    const uint32_t tb_B = graph.get_id("tb_B");
    const uint32_t tb_C = graph.get_id("tb_C");

    EXPECT_EQ(r.phasing_nodes.count(tb_A), 1u);
    EXPECT_EQ(r.phasing_nodes.count(tb_B), 1u);
    EXPECT_EQ(r.phasing_nodes.count(tb_C), 1u);

    EXPECT_EQ(r.alt_map.count(tb_A), 0u);
    EXPECT_EQ(r.alt_map.count(tb_B), 0u);
    EXPECT_EQ(r.alt_map.count(tb_C), 0u);
}

TEST_F(UnifiedBubbleTest, NestedOuterBubble) {
    const uint32_t nb_E         = graph.get_id("nb_E");
    const uint32_t nb_inner_src = graph.get_id("nb_inner_src");
    const uint32_t nb_inner_snk = graph.get_id("nb_inner_snk");

    EXPECT_EQ(r.phasing_nodes.count(nb_E),         1u);
    EXPECT_EQ(r.phasing_nodes.count(nb_inner_src), 1u);
    EXPECT_EQ(r.phasing_nodes.count(nb_inner_snk), 1u);

    EXPECT_EQ(r.alt_map.at(nb_E),         nb_inner_src);
    EXPECT_EQ(r.alt_map.at(nb_inner_src), nb_E);
    EXPECT_EQ(r.alt_map.at(nb_inner_snk), nb_E);
}

TEST_F(UnifiedBubbleTest, NestedInnerBubble) {
    const uint32_t ni_A = graph.get_id("ni_A");
    const uint32_t ni_B = graph.get_id("ni_B");

    EXPECT_EQ(r.phasing_nodes.count(ni_A), 1u);
    EXPECT_EQ(r.phasing_nodes.count(ni_B), 1u);

    EXPECT_EQ(r.alt_map.at(ni_A), ni_B);
    EXPECT_EQ(r.alt_map.at(ni_B), ni_A);
}

TEST_F(UnifiedBubbleTest, StarSequenceBubbleStillRecorded) {
    const uint32_t star_A = graph.get_id("star_A");
    const uint32_t star_B = graph.get_id("star_B");

    EXPECT_EQ(r.phasing_nodes.count(star_A), 1u);
    EXPECT_EQ(r.phasing_nodes.count(star_B), 1u);
    EXPECT_EQ(r.alt_map.at(star_A), star_B);
    EXPECT_EQ(r.alt_map.at(star_B), star_A);
}

TEST_F(UnifiedBubbleTest, LinearChainNotPhased) {
    for (const char* name : {"lc_a", "lc_b", "lc_c", "lc_d"}) {
        EXPECT_EQ(r.phasing_nodes.count(graph.get_id(name)), 0u)
            << name << " should not be a phasing node";
    }
}

TEST_F(UnifiedBubbleTest, StructuralNodesNotPhased) {
    for (const char* name : {
            "bb_start", "bb_end",
            "snp_src",  "snp_snk",
            "mb_src",   "mb_snk",
            "tb_src",   "tb_snk",
            "nb_src",   "nb_snk",
            "sim_src",  "sim_snk",
            "dis_src",  "dis_snk",
            "star_src", "star_snk"}) {
        EXPECT_EQ(r.phasing_nodes.count(graph.get_id(name)), 0u)
            << name << " should not be a phasing node";
    }
}

TEST_F(UnifiedBubbleTest, PhasingNodeCount) {
    EXPECT_EQ(r.phasing_nodes.size(), 20u);
}

TEST_F(UnifiedTest, MinHashFiltersDissimilarBranches) {
    auto rfilt = det.find_unlabeled_alts(graph, 16, 1000, 0.5);

    const uint32_t sim_A = graph.get_id("sim_A");
    const uint32_t sim_B = graph.get_id("sim_B");
    const uint32_t dis_A = graph.get_id("dis_A");
    const uint32_t dis_B = graph.get_id("dis_B");

    EXPECT_EQ(rfilt.phasing_nodes.count(sim_A), 1u);
    EXPECT_EQ(rfilt.phasing_nodes.count(sim_B), 1u);

    EXPECT_EQ(rfilt.phasing_nodes.count(dis_A), 0u);
    EXPECT_EQ(rfilt.phasing_nodes.count(dis_B), 0u);
}

TEST_F(UnifiedTest, ShastaPRStylePair) {
    auto rs = det.get_alts_from_shasta_names(graph);

    const uint32_t pr0 = graph.get_id("PR.000001.0");
    const uint32_t pr1 = graph.get_id("PR.000001.1");

    EXPECT_EQ(rs.phasing_nodes.count(pr0), 1u);
    EXPECT_EQ(rs.phasing_nodes.count(pr1), 1u);
    EXPECT_EQ(rs.alt_map.at(pr0), pr1);
    EXPECT_EQ(rs.alt_map.at(pr1), pr0);
}

TEST_F(UnifiedTest, ShastaDotSuffixPair) {
    auto rs = det.get_alts_from_shasta_names(graph);

    const uint32_t c0 = graph.get_id("contig1.0");
    const uint32_t c1 = graph.get_id("contig1.1");

    EXPECT_EQ(rs.phasing_nodes.count(c0), 1u);
    EXPECT_EQ(rs.phasing_nodes.count(c1), 1u);
    EXPECT_EQ(rs.alt_map.at(c0), c1);
    EXPECT_EQ(rs.alt_map.at(c1), c0);
}

TEST_F(UnifiedTest, ShastaGroupOfThreeSkipped) {
    auto rs = det.get_alts_from_shasta_names(graph);

    EXPECT_EQ(rs.phasing_nodes.count(graph.get_id("base.0")), 0u);
    EXPECT_EQ(rs.phasing_nodes.count(graph.get_id("base.1")), 0u);
    EXPECT_EQ(rs.phasing_nodes.count(graph.get_id("base.2")), 0u);
}

TEST_F(UnifiedTest, ShastaNoMatchingAlt) {
    auto rs = det.get_alts_from_shasta_names(graph);
    EXPECT_EQ(rs.phasing_nodes.count(graph.get_id("uniquenode")), 0u);
}

TEST_F(UnifiedTest, ShastaClassicNamesIgnoredUnderDotSuffixStyle) {
    auto rs = det.get_alts_from_shasta_names(graph);

    EXPECT_EQ(rs.phasing_nodes.count(graph.get_id("000001F")),     0u);
    EXPECT_EQ(rs.phasing_nodes.count(graph.get_id("000001F_003")), 0u);
}

TEST(BubbleDetectorEdgeCases, EmptyGraph) {
    Graph g;
    BubbleDetector det;
    auto rs = det.find_unlabeled_alts(g);
    EXPECT_TRUE(rs.phasing_nodes.empty());
    EXPECT_TRUE(rs.alt_map.empty());
}

TEST_F(UnifiedTest, ContactFilterPreservesBubblePairsRemovesOthers) {
    auto rb = det.find_unlabeled_alts(graph, 16, 1000, 0.0);

    ContactMatrix cm;
    cm.load_contacts(std::string(TEST_DATA_DIR) + "/test_unified.bam", graph);
    cm.build_csr(graph.get_num_nodes());

    const uint32_t lc_a     = graph.get_id("lc_a");
    const uint32_t lc_b     = graph.get_id("lc_b");
    const uint32_t bb_start = graph.get_id("bb_start");
    const uint32_t snp_A    = graph.get_id("snp_A");
    const uint32_t snp_B    = graph.get_id("snp_B");
    const uint32_t mb_A     = graph.get_id("mb_A");

    EXPECT_EQ(cm.get_contact(snp_A, snp_B),    3u);
    EXPECT_EQ(cm.get_contact(snp_A, mb_A),     1u);
    EXPECT_EQ(cm.get_contact(snp_A, bb_start), 1u);
    EXPECT_EQ(cm.get_contact(lc_a,  lc_b),     1u);

    cm.filter_to_phasing_nodes(rb.phasing_nodes, graph.get_num_nodes());

    EXPECT_EQ(cm.get_contact(snp_A, snp_B), 3u);
    EXPECT_EQ(cm.get_contact(snp_A, mb_A),  1u);

    EXPECT_EQ(cm.get_contact(snp_A, bb_start), 0u);
    EXPECT_EQ(cm.get_contact(lc_a,  lc_b),     0u);
}

TEST(ContactMatrixFilterTest, DropsNonPhasingContacts) {
    ContactMatrix cm;
    cm.add_contact(0, 1, 7);  
    cm.add_contact(0, 2, 4);  
    cm.add_contact(1, 2, 2); 
    cm.add_contact(0, 3, 9);  
    cm.add_contact(2, 4, 6);  
    cm.add_contact(3, 4, 5);  
    cm.add_contact(4, 5, 1);  
    cm.build_csr(6);

    EXPECT_EQ(cm.get_contact(0, 1), 7u);
    EXPECT_EQ(cm.get_contact(0, 3), 9u);
    EXPECT_EQ(cm.get_contact(3, 4), 5u);

    std::unordered_set<uint32_t> phasing = {0, 1, 2};
    cm.filter_to_phasing_nodes(phasing, 6);

    EXPECT_EQ(cm.get_contact(0, 1), 7u);
    EXPECT_EQ(cm.get_contact(0, 2), 4u);
    EXPECT_EQ(cm.get_contact(1, 2), 2u);

    EXPECT_EQ(cm.get_contact(0, 3), 0u);
    EXPECT_EQ(cm.get_contact(2, 4), 0u);
    EXPECT_EQ(cm.get_contact(3, 4), 0u);
    EXPECT_EQ(cm.get_contact(4, 5), 0u);

    EXPECT_EQ(cm.row_begin(3), cm.row_end(3));
    EXPECT_EQ(cm.row_begin(4), cm.row_end(4));
}

TEST_F(UnifiedTest, TrioConstraintsLockDominantNodes) {
    auto rb = det.find_unlabeled_alts(graph, 16, 1000, 0.0);

    const uint32_t snp_A  = graph.get_id("snp_A");   
    const uint32_t snp_B  = graph.get_id("snp_B");   
    const uint32_t mb_A   = graph.get_id("mb_A");    
    const uint32_t mb_B   = graph.get_id("mb_B");    

    std::vector<TrioScores> scores(graph.get_num_nodes());
    scores[snp_A] = {100, 5};    
    scores[snp_B] = {5, 100};    
    scores[mb_A]  = {50, 50};    
    scores[mb_B]  = {0, 0};      

    size_t locked = apply_trio_constraints(graph, rb.phasing_nodes, scores, PhaserConfig{});
    EXPECT_EQ(locked, 2u);

    EXPECT_TRUE (graph.is_phase_locked(snp_A));
    EXPECT_EQ   (graph.get_phase      (snp_A), 0);
    EXPECT_TRUE (graph.is_phase_locked(snp_B));
    EXPECT_EQ   (graph.get_phase      (snp_B), 1);

    EXPECT_FALSE(graph.is_phase_locked(mb_A));
    EXPECT_EQ   (graph.get_phase      (mb_A), -1);
    EXPECT_FALSE(graph.is_phase_locked(mb_B));
    EXPECT_EQ   (graph.get_phase      (mb_B), -1);
}

