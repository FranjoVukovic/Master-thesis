#pragma once

#include <memory>
#include <thread_pool/thread_pool.hpp>

#include "Chainer.hpp"

class Graph;
class ContactMatrix;

class HamiltonianChainer {
public:
    explicit HamiltonianChainer(std::shared_ptr<thread_pool::ThreadPool> pool);

    ChainResult generate_chain_paths(const Graph& graph,
                                     const BubbleResult& bubbles,
                                     const ContactMatrix& contacts) const;

private:
    std::shared_ptr<thread_pool::ThreadPool> pool_;
};
