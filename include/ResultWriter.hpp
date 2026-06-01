#pragma once

#include <memory>
#include <string>
#include <thread_pool/thread_pool.hpp>

#include "Chainer.hpp"
#include "UnzippedGraph.hpp"

class Graph;

struct OutputConfig {
    std::string out_dir = ".";
    std::string prefix  = "out";
};

class ResultWriter {
public:
    explicit ResultWriter(std::shared_ptr<thread_pool::ThreadPool> pool = nullptr);

    void write_chained_gfa (const Graph& graph,
                            const ChainResult& chains,
                            const OutputConfig& cfg) const;

    void write_unzipped_gfa(const Graph& graph,
                            const UnzippedGraph& uz,
                            const OutputConfig& cfg) const;

    // One concatenated sequence per chain (haplotype path), mirroring how
    // GFAse emits one contig per walked path — not one record per node.
    void write_fastas      (const Graph& graph,
                            const ChainResult& chains,
                            const OutputConfig& cfg) const;

private:
    std::shared_ptr<thread_pool::ThreadPool> pool_;
};
