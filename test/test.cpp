#include <gtest/gtest.h>

#include "BubbleDetector.hpp"
#include "ContactMatrix.hpp"
#include "Graph.hpp"
#include "Phaser.hpp"
#include "TrioBinner.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread_pool/thread_pool.hpp>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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

// ============================================================
// Phaser test suite (T01..T14)
// ============================================================

namespace {
const std::string kGfaPath = std::string(TEST_DATA_DIR) + "/test_unified.gfa";
const std::string kBamPath = std::string(TEST_DATA_DIR) + "/test_unified.bam";
const std::string kPatYak  = std::string(TEST_DATA_DIR) + "/test_unified_pat.yak";
const std::string kMatYak  = std::string(TEST_DATA_DIR) + "/test_unified_mat.yak";

bool yak_files_exist() {
    std::ifstream a(kPatYak), b(kMatYak);
    return a.good() && b.good();
}

struct Pipeline {
    Graph         graph;
    ContactMatrix contacts;
    BubbleResult  bubbles;
};

void setup_pipeline(Pipeline& p) {
    p.graph.load_from_gfa(kGfaPath);
    p.contacts.load_contacts(kBamPath, p.graph, 1);
    p.contacts.build_csr(p.graph.get_num_nodes(), nullptr);
    BubbleDetector det;
    p.bubbles = det.find_unlabeled_alts(p.graph, 16, 1000, 0.2, nullptr);
    p.contacts.filter_to_phasing_nodes(p.bubbles.phasing_nodes,
                                        p.graph.get_num_nodes());
}
} // namespace

// ---------- T01 Graph loading ----------
TEST(PhaserT01_Graph, LoadingAndAccessors) {
    Graph g;
    ASSERT_NO_THROW(g.load_from_gfa(kGfaPath));
    EXPECT_EQ(g.get_num_nodes(), 50u);

    EXPECT_NO_THROW((void) g.get_id("snp_A"));
    EXPECT_EQ(g.get_name(g.get_id("snp_A")), "snp_A");
    EXPECT_THROW((void) g.get_id("nonexistent_node_xyz"), std::out_of_range);

    const uint32_t id = g.get_id("snp_A");
    EXPECT_EQ(g.get_phase(id), -1);
    EXPECT_FALSE(g.is_locked(id));
}

// ---------- T02 Contact matrix loading + CSR accessors ----------
TEST(PhaserT02_Contact, LoadAndCsrAccessors) {
    Graph g;
    g.load_from_gfa(kGfaPath);
    ContactMatrix cm;
    cm.load_contacts(kBamPath, g, 1);
    cm.build_csr(g.get_num_nodes(), nullptr);

    const uint32_t snp_A = g.get_id("snp_A");
    const uint32_t snp_B = g.get_id("snp_B");
    const uint32_t mb_A  = g.get_id("mb_A");
    const uint32_t mb_B  = g.get_id("mb_B");
    const uint32_t lc_a  = g.get_id("lc_a");
    const uint32_t lc_b  = g.get_id("lc_b");

    EXPECT_EQ(cm.get_contact(snp_A, snp_B), 3u);
    EXPECT_EQ(cm.get_contact(snp_B, snp_A), 3u);
    EXPECT_EQ(cm.get_contact(snp_A, mb_A),  1u);
    EXPECT_EQ(cm.get_contact(mb_A,  snp_A), 1u);
    EXPECT_EQ(cm.get_contact(lc_a,  lc_b),  1u);
    EXPECT_EQ(cm.get_contact(snp_A, mb_B),  0u);

    bool found_snpB = false, found_mbA = false;
    bool wgt_snpB_ok = false, wgt_mbA_ok = false;
    for (size_t pos = cm.row_begin(snp_A); pos != cm.row_end(snp_A); ++pos) {
        const uint32_t col = cm.col_at(pos);
        const uint32_t val = cm.value_at(pos);
        if (col == snp_B) { found_snpB = true; wgt_snpB_ok = (val == 3); }
        if (col == mb_A)  { found_mbA  = true; wgt_mbA_ok  = (val == 1); }
    }
    EXPECT_TRUE(found_snpB);
    EXPECT_TRUE(wgt_snpB_ok);
    EXPECT_TRUE(found_mbA);
    EXPECT_TRUE(wgt_mbA_ok);
}

// ---------- T03 filter_to_phasing_nodes ----------
TEST(PhaserT03_Contact, FilterToPhasingNodes) {
    Graph g;
    g.load_from_gfa(kGfaPath);
    ContactMatrix cm;
    cm.load_contacts(kBamPath, g, 1);
    cm.build_csr(g.get_num_nodes(), nullptr);

    BubbleDetector det;
    BubbleResult br = det.find_unlabeled_alts(g, 16, 1000, 0.2, nullptr);
    cm.filter_to_phasing_nodes(br.phasing_nodes, g.get_num_nodes());

    const uint32_t lc_a     = g.get_id("lc_a");
    const uint32_t lc_b     = g.get_id("lc_b");
    const uint32_t snp_A    = g.get_id("snp_A");
    const uint32_t snp_B    = g.get_id("snp_B");
    const uint32_t mb_A     = g.get_id("mb_A");
    const uint32_t bb_start = g.get_id("bb_start");

    EXPECT_EQ(cm.get_contact(lc_a, lc_b),     0u);
    EXPECT_EQ(cm.get_contact(snp_A, bb_start), 0u);
    EXPECT_EQ(cm.get_contact(snp_A, snp_B),    3u);
    EXPECT_EQ(cm.get_contact(snp_A, mb_A),     1u);
}

// ---------- T04 Superbubble mode ----------
TEST(PhaserT04_BubbleDet, SuperbubbleMode) {
    Graph g;
    g.load_from_gfa(kGfaPath);
    BubbleDetector det;
    BubbleResult br = det.find_unlabeled_alts(g, 16, 1000, 0.2, nullptr);

    const uint32_t snp_A   = g.get_id("snp_A");
    const uint32_t snp_B   = g.get_id("snp_B");
    const uint32_t ni_A    = g.get_id("ni_A");
    const uint32_t ni_B    = g.get_id("ni_B");
    const uint32_t sim_A   = g.get_id("sim_A");
    const uint32_t sim_B   = g.get_id("sim_B");
    const uint32_t star_A  = g.get_id("star_A");
    const uint32_t star_B  = g.get_id("star_B");
    const uint32_t snp_src = g.get_id("snp_src");
    const uint32_t snp_snk = g.get_id("snp_snk");
    const uint32_t uniqn   = g.get_id("uniquenode");

    EXPECT_EQ(br.phasing_nodes.count(snp_A), 1u);
    EXPECT_EQ(br.phasing_nodes.count(snp_B), 1u);
    EXPECT_EQ(br.alt_map.at(snp_A), snp_B);
    EXPECT_EQ(br.alt_map.at(snp_B), snp_A);

    EXPECT_EQ(br.phasing_nodes.count(ni_A), 1u);
    EXPECT_EQ(br.phasing_nodes.count(ni_B), 1u);
    EXPECT_EQ(br.alt_map.at(ni_A), ni_B);

    EXPECT_EQ(br.phasing_nodes.count(snp_src), 0u);
    EXPECT_EQ(br.phasing_nodes.count(snp_snk), 0u);
    EXPECT_EQ(br.phasing_nodes.count(uniqn),   0u);

    EXPECT_EQ(br.alt_map.count(sim_A), 1u);
    EXPECT_EQ(br.alt_map.at(sim_A),    sim_B);

    // Star bubble: both seqs empty → detector skips Jaccard filter
    // (find_unlabeled_alts only filters when BOTH seqs are non-empty),
    // so the pair is still recorded based on topology alone.
    EXPECT_EQ(br.alt_map.count(star_A), 1u);
    EXPECT_EQ(br.alt_map.at(star_A),    star_B);

    // dis bubble: both seqs non-empty (AAAA vs CCCC) → Jaccard ≈ 0 < 0.2 → filtered out.
    const uint32_t dis_A = g.get_id("dis_A");
    const uint32_t dis_B = g.get_id("dis_B");
    EXPECT_EQ(br.alt_map.count(dis_A), 0u);
    EXPECT_EQ(br.alt_map.count(dis_B), 0u);
}

// ---------- T05 Shasta name mode ----------
TEST(PhaserT05_BubbleDet, ShastaNameMode) {
    Graph g;
    g.load_from_gfa(kGfaPath);
    BubbleDetector det;
    BubbleResult br = det.get_alts_from_shasta_names(g);

    const uint32_t pr0 = g.get_id("PR.000001.0");
    const uint32_t pr1 = g.get_id("PR.000001.1");
    const uint32_t c0  = g.get_id("contig1.0");
    const uint32_t c1  = g.get_id("contig1.1");
    const uint32_t f0  = g.get_id("000001F");
    const uint32_t f1  = g.get_id("000001F_003");
    const uint32_t b0  = g.get_id("base.0");
    const uint32_t b1  = g.get_id("base.1");
    const uint32_t b2  = g.get_id("base.2");
    const uint32_t uniqn = g.get_id("uniquenode");

    auto pair_ok = [&](uint32_t a, uint32_t b) {
        return br.alt_map.count(a) == 1 && br.alt_map.count(b) == 1 &&
               br.alt_map.at(a) == b && br.alt_map.at(b) == a;
    };
    EXPECT_TRUE(pair_ok(pr0, pr1));
    EXPECT_TRUE(pair_ok(c0,  c1));

    // detect_style picks one assembler-naming style per graph; with PR/dot-suffix
    // patterns dominant in the unified GFA, classic-Shasta names like 000001F /
    // 000001F_003 are intentionally NOT paired.
    EXPECT_FALSE(pair_ok(f0, f1));
    EXPECT_EQ(br.phasing_nodes.count(f0), 0u);
    EXPECT_EQ(br.phasing_nodes.count(f1), 0u);

    EXPECT_EQ(br.phasing_nodes.count(b0), 0u);
    EXPECT_EQ(br.phasing_nodes.count(b1), 0u);
    EXPECT_EQ(br.phasing_nodes.count(b2), 0u);
    EXPECT_EQ(br.phasing_nodes.count(uniqn), 0u);
}

// ---------- T06 apply_trio_constraints with real yak ----------
TEST(PhaserT06_Trio, ApplyConstraints) {
    if (!yak_files_exist()) GTEST_SKIP() << "yak files unavailable";

    Graph g;
    g.load_from_gfa(kGfaPath);
    auto pool = std::make_shared<thread_pool::ThreadPool>(2);

    std::unique_ptr<TrioBinner> binner;
    try { binner = std::make_unique<TrioBinner>(kPatYak, kMatYak); }
    catch (const std::exception&) { GTEST_SKIP() << "yak files unavailable"; }

    const std::vector<TrioScores> scores = binner->compute_scores(g, pool);

    const uint32_t dis_A = g.get_id("dis_A");
    const uint32_t dis_B = g.get_id("dis_B");
    const uint32_t snp_A = g.get_id("snp_A");

    const uint64_t dA_total = uint64_t{scores[dis_A].pat_count} + scores[dis_A].mat_count;
    const uint64_t dB_total = uint64_t{scores[dis_B].pat_count} + scores[dis_B].mat_count;
    if (dA_total == 0 || dB_total == 0) {
        GTEST_SKIP() << "yak DBs produce zero counts on this graph (k-mer mismatch)";
    }

    std::unordered_set<uint32_t> phasing_nodes = {dis_A, dis_B, snp_A, g.get_id("snp_B")};

    PhaserConfig cfg;
    cfg.lock_ratio      = 0.90;
    cfg.min_total_kmers = 1;
    apply_trio_constraints(g, phasing_nodes, scores, cfg);

    EXPECT_TRUE(g.is_locked(dis_A));
    EXPECT_EQ  (g.get_phase(dis_A), 0);
    EXPECT_TRUE(g.is_locked(dis_B));
    EXPECT_EQ  (g.get_phase(dis_B), 1);
    EXPECT_FALSE(g.is_locked(snp_A));
    EXPECT_EQ  (g.get_phase(snp_A), -1);
}

// ---------- T07 monte_carlo_phase correctness on snp bubble ----------
TEST(PhaserT07_MC, SnpBubbleCorrectness) {
    Pipeline p;
    setup_pipeline(p);

    PhaserConfig cfg;
    cfg.core_iterations = 200;
    cfg.n_rounds        = 2;
    cfg.rng_seed        = 42;

    auto pool = std::make_shared<thread_pool::ThreadPool>(4);
    monte_carlo_phase(p.graph, p.contacts, p.bubbles, cfg, pool);

    const uint32_t snp_A = p.graph.get_id("snp_A");
    const uint32_t snp_B = p.graph.get_id("snp_B");
    const uint32_t mb_A  = p.graph.get_id("mb_A");
    const uint32_t mb_B  = p.graph.get_id("mb_B");

    EXPECT_NE(p.graph.get_phase(snp_A), -1);
    EXPECT_NE(p.graph.get_phase(snp_B), -1);
    EXPECT_EQ(p.graph.get_phase(snp_A) + p.graph.get_phase(snp_B), 1);
    EXPECT_EQ(p.graph.get_phase(snp_A), p.graph.get_phase(mb_A));
    EXPECT_EQ(p.graph.get_phase(snp_B), p.graph.get_phase(mb_B));

    for (uint32_t id : p.bubbles.phasing_nodes) {
        if (!p.bubbles.alt_map.count(id)) continue;
        const int8_t ph = p.graph.get_phase(id);
        EXPECT_TRUE(ph == 0 || ph == 1)
            << "node " << p.graph.get_name(id) << " has phase " << static_cast<int>(ph);
    }
}

// ---------- T08 MC: trio-locked nodes never overwritten ----------
TEST(PhaserT08_MC, TrioLockedPreserved) {
    if (!yak_files_exist()) GTEST_SKIP() << "yak files unavailable";

    Pipeline p;
    setup_pipeline(p);
    auto pool = std::make_shared<thread_pool::ThreadPool>(4);

    std::unique_ptr<TrioBinner> binner;
    try { binner = std::make_unique<TrioBinner>(kPatYak, kMatYak); }
    catch (const std::exception&) { GTEST_SKIP() << "yak files unavailable"; }

    const auto scores = binner->compute_scores(p.graph, pool);

    const uint32_t dis_A = p.graph.get_id("dis_A");
    const uint32_t dis_B = p.graph.get_id("dis_B");

    const uint64_t dA_total = uint64_t{scores[dis_A].pat_count} + scores[dis_A].mat_count;
    const uint64_t dB_total = uint64_t{scores[dis_B].pat_count} + scores[dis_B].mat_count;
    if (dA_total == 0 || dB_total == 0) {
        GTEST_SKIP() << "yak DBs produce zero counts on this graph (k-mer mismatch)";
    }

    PhaserConfig cfg;
    cfg.lock_ratio      = 0.90;
    cfg.min_total_kmers = 1;
    apply_trio_constraints(p.graph, p.bubbles.phasing_nodes, scores, cfg);

    monte_carlo_phase(p.graph, p.contacts, p.bubbles, cfg, pool);

    EXPECT_EQ (p.graph.get_phase(dis_A), 0);
    EXPECT_EQ (p.graph.get_phase(dis_B), 1);
    EXPECT_TRUE(p.graph.is_locked(dis_A));
    EXPECT_TRUE(p.graph.is_locked(dis_B));
}

// ---------- T09 MC: empty bubble list returns immediately ----------
TEST(PhaserT09_MC, EmptyBubblesNoChange) {
    Graph g;
    g.load_from_gfa(kGfaPath);
    ContactMatrix cm;
    cm.load_contacts(kBamPath, g, 1);
    cm.build_csr(g.get_num_nodes(), nullptr);

    BubbleResult empty;
    PhaserConfig cfg;
    auto pool = std::make_shared<thread_pool::ThreadPool>(2);

    EXPECT_NO_THROW(monte_carlo_phase(g, cm, empty, cfg, pool));

    for (uint32_t i = 0; i < g.get_num_nodes(); ++i) {
        EXPECT_EQ(g.get_phase(i), -1);
    }
}

// ---------- T10 MC: all bubbles trio-locked → no free pairs ----------
TEST(PhaserT10_MC, AllLockedNoChange) {
    Pipeline p;
    setup_pipeline(p);
    auto pool = std::make_shared<thread_pool::ThreadPool>(2);

    std::unordered_map<uint32_t, int8_t> snapshot;
    int flip = 0;
    for (const auto& kv : p.bubbles.alt_map) {
        const uint32_t a = kv.first, b = kv.second;
        if (a >= b) continue;
        const int8_t pa = (flip & 1) ? int8_t{1} : int8_t{0};
        const int8_t pb = static_cast<int8_t>(1 - pa);
        p.graph.lock_phase(a, pa);
        p.graph.lock_phase(b, pb);
        snapshot[a] = pa;
        snapshot[b] = pb;
        ++flip;
    }

    EXPECT_NO_THROW(monte_carlo_phase(p.graph, p.contacts, p.bubbles,
                                       PhaserConfig{}, pool));
    for (const auto& kv : snapshot) {
        EXPECT_EQ (p.graph.get_phase(kv.first), kv.second);
        EXPECT_TRUE(p.graph.is_locked(kv.first));
    }
}

// ---------- T11 MC: reproducibility ----------
TEST(PhaserT11_MC, Reproducibility) {
    auto run_once = [](std::unordered_map<std::string, int8_t>& out) {
        Pipeline p;
        setup_pipeline(p);
        PhaserConfig cfg; cfg.rng_seed = 42;
        auto pool = std::make_shared<thread_pool::ThreadPool>(4);
        monte_carlo_phase(p.graph, p.contacts, p.bubbles, cfg, pool);
        for (uint32_t id : p.bubbles.phasing_nodes) {
            if (!p.bubbles.alt_map.count(id)) continue;
            out[p.graph.get_name(id)] = p.graph.get_phase(id);
        }
    };

    std::unordered_map<std::string, int8_t> a, b;
    run_once(a);
    run_once(b);

    ASSERT_EQ(a.size(), b.size());
    for (const auto& kv : a) {
        ASSERT_EQ(b.count(kv.first), 1u) << kv.first;
        EXPECT_EQ(b.at(kv.first), kv.second) << kv.first;
    }
}

// ---------- T12 MC: thread count scales sample count ----------
TEST(PhaserT12_MC, ThreadScalingPreservesInvariants) {
    auto run_with = [](size_t n_threads,
                       int8_t& snp_A_ph, int8_t& snp_B_ph,
                       bool& all_valid) {
        Pipeline p;
        setup_pipeline(p);
        PhaserConfig cfg; cfg.rng_seed = 42;
        auto pool = std::make_shared<thread_pool::ThreadPool>(n_threads);
        monte_carlo_phase(p.graph, p.contacts, p.bubbles, cfg, pool);
        snp_A_ph = p.graph.get_phase(p.graph.get_id("snp_A"));
        snp_B_ph = p.graph.get_phase(p.graph.get_id("snp_B"));
        all_valid = true;
        for (uint32_t id : p.bubbles.phasing_nodes) {
            if (!p.bubbles.alt_map.count(id)) continue;
            const int8_t ph = p.graph.get_phase(id);
            if (ph != 0 && ph != 1) { all_valid = false; break; }
        }
    };

    int8_t a2 = 0, b2 = 0, a8 = 0, b8 = 0;
    bool valid_2 = false, valid_8 = false;
    run_with(2, a2, b2, valid_2);
    run_with(8, a8, b8, valid_8);

    EXPECT_TRUE(valid_2);
    EXPECT_TRUE(valid_8);
    EXPECT_EQ(a2 + b2, 1);
    EXPECT_EQ(a8 + b8, 1);
}

// ---------- T13 ΔS sign check (in-memory unit) ----------
TEST(PhaserT13_MC, UnitSimpleBubble) {
    Graph g;
    g.load_from_gfa_string(
        "S\tA\tACGT\nS\tB\tACGT\nS\tC\tACGT\nS\tD\tACGT\n"
        "L\tA\t+\tB\t+\t0M\nL\tA\t+\tC\t+\t0M\n"
        "L\tB\t+\tD\t+\t0M\nL\tC\t+\tD\t+\t0M\n");

    const uint32_t idB = g.get_id("B");
    const uint32_t idC = g.get_id("C");

    ContactMatrix cm;
    cm.add_contact(idB, idC, 10);
    cm.build_csr(g.get_num_nodes(), nullptr);

    BubbleResult br;
    br.phasing_nodes.insert(idB);
    br.phasing_nodes.insert(idC);
    br.alt_map[idB] = idC;
    br.alt_map[idC] = idB;
    cm.filter_to_phasing_nodes(br.phasing_nodes, g.get_num_nodes());

    PhaserConfig cfg;
    cfg.rng_seed        = 0;
    cfg.n_rounds        = 1;
    cfg.core_iterations = 50;
    auto pool = std::make_shared<thread_pool::ThreadPool>(2);
    monte_carlo_phase(g, cm, br, cfg, pool);

    EXPECT_NE(g.get_phase(idB), -1);
    EXPECT_NE(g.get_phase(idC), -1);
    EXPECT_EQ(g.get_phase(idB) + g.get_phase(idC), 1);
}

// ---------- T14 CSV output ----------
TEST(PhaserT14_CSV, OutputContent) {
    Pipeline p;
    setup_pipeline(p);
    PhaserConfig cfg; cfg.rng_seed = 42;
    auto pool = std::make_shared<thread_pool::ThreadPool>(4);
    monte_carlo_phase(p.graph, p.contacts, p.bubbles, cfg, pool);

    const std::string csv_path = std::string(TEST_DATA_DIR) + "/phase_labels_test.csv";
    {
        std::ofstream out(csv_path);
        out << "contig_id,phase\n";
        for (uint32_t id = 0; id < p.graph.get_num_nodes(); ++id) {
            const int ph = p.graph.get_phase(id);
            if (ph == -1) continue;
            out << p.graph.get_name(id) << "," << ph << "\n";
        }
    }

    std::ifstream in(csv_path);
    ASSERT_TRUE(in.good()) << "cannot open " << csv_path;
    std::string header;
    std::getline(in, header);
    EXPECT_EQ(header, "contig_id,phase");

    std::unordered_map<std::string, int> rows;
    std::string line;
    bool only_01 = true;
    bool no_dup  = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const size_t comma = line.find(',');
        ASSERT_NE(comma, std::string::npos) << "bad row: " << line;
        const std::string name = line.substr(0, comma);
        const int ph = std::stoi(line.substr(comma + 1));
        if (ph != 0 && ph != 1) only_01 = false;
        if (rows.count(name))   no_dup  = false;
        rows[name] = ph;
    }
    in.close();
    std::remove(csv_path.c_str());

    EXPECT_TRUE(only_01);
    EXPECT_TRUE(no_dup);
    ASSERT_EQ(rows.count("snp_A"), 1u);
    ASSERT_EQ(rows.count("snp_B"), 1u);
    EXPECT_EQ(rows["snp_A"] + rows["snp_B"], 1);
    EXPECT_EQ(rows.count("lc_a"),       0u);
    EXPECT_EQ(rows.count("snp_src"),    0u);
    EXPECT_EQ(rows.count("uniquenode"), 0u);
}

