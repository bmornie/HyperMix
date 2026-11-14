#include "generator.hpp"

#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <numeric>

#include "random.hpp"
#include "hypergraph.hpp"



// ============================================================================
// Hyperedge formation mechanisms
// ============================================================================

int random_edge(Hypergraph& H, size_t edge_size, double p_new, int max_attempts)
{
    std::vector<ID> new_edge {};
    new_edge.reserve(edge_size);

    ID N {static_cast<ID>(H.node_count())};
    for (size_t k {}; k < edge_size; ++k) {
        if (N < edge_size || random_uniform() < p_new) {
            new_edge.push_back(N++);
        }
    }

    edge_size -= new_edge.size();

    if (edge_size == 0) {
        H.add_edge(new_edge);
        return 1;
    }

    for (int i {}; i < max_attempts; ++i)
    {
        std::vector<ID> nodes {rand_int((ID)0, static_cast<ID>(H.node_count())-1, edge_size)};

        nodes.insert(nodes.end(), new_edge.begin(), new_edge.end());

        if (H.has_edge(nodes))
            continue;

        H.add_edge(nodes);
        return 1;
    }
    return 0;
}

int PA(Hypergraph& H, size_t edge_size, double p_new, int max_attempts)
{
    if (H.edge_count() == 0) { return 0; }

    ID N {static_cast<ID>(H.node_count())};
    ID seed_node {rand_int((ID)0, N-1)};
    if (random_uniform() < p_new) { seed_node = N; }

    const auto nodes {H.nodes()};
    std::vector<double> degrees (nodes.size());
    size_t larger_than_zero {};
    for (size_t i {}; i < degrees.size(); ++i)
    {
        auto d {H.degree(nodes[i])};
        if (nodes[i] != seed_node) {
            degrees[i] = static_cast<double>(d);
            if (d > 0) { ++larger_than_zero; }
        }
        else {
            degrees[i] = 0;     // ensure seed_node will not be included twice
        }
    }

    if (edge_size-1 > larger_than_zero)
        return 0;

    for (int i {}; i < max_attempts; ++i)
    {
        auto new_edge = weighted_sample_without_replacement(nodes, degrees, edge_size-1);
        new_edge.push_back(seed_node);

        if (H.has_edge(new_edge))
            continue;

        H.add_edge(new_edge);
        return 1;
    }
    return 0;
}

int mutation(Hypergraph &H, double p_new, int max_attempts)
{
    if (H.edge_count() == 0) { return 0; }

    for (int i {}; i < max_attempts; ++i) {
        auto new_edge {H.get_edge(get_sample(H.edges()))};

        // remove a node
        pop_random(new_edge);

        // add a new node

        ID N {static_cast<ID>(H.node_count())};
        if (new_edge.size() == N || random_uniform() < p_new) {
            new_edge.push_back(N);
            H.add_edge(new_edge);
            return 1;
        }

        ID new_node {};
        while (true)
        {
            new_node = rand_int((ID)0, N-1);
            if (std::find(new_edge.begin(), new_edge.end(), new_node) == new_edge.end())
                break;
        }

        new_edge.push_back(new_node);

        if (H.has_edge(new_edge))
            continue;
        
        H.add_edge(new_edge);
        return 1;
    }
    return 0;
}

int nesting(Hypergraph& H, int max_attempts)
{
    std::vector<ID> large_edges {};
    for (auto e : H.edges()) {
        if (H.size(e) >= 3) {
            large_edges.push_back(e);
        }
    }

    if (large_edges.size() == 0) { return 0; }

    for (int i {}; i < max_attempts; ++i)
    {
        auto old_edge {H.get_edge(get_sample(large_edges))};
        size_t new_size {rand_int((size_t)2, old_edge.size() - 1)};
        auto new_edge {sample_without_replacement(old_edge, new_size)};

        if (H.has_edge(new_edge))
            continue;
        
        H.add_edge(new_edge);
        return 1;
    }
    return 0;
}

int merging(Hypergraph& H, int max_attempts)
{
    if (H.edge_count() < 2) { return 0; }

    const auto edges {H.edges()};

    for (int i {}; i < max_attempts; ++i) {
        auto merge_ids {sample_without_replacement(edges, 2)};

        auto new_edge {H.get_edge(merge_ids[0])};

        for (auto n : H.get_edge(merge_ids[1])){
            if (std::find(new_edge.begin(), new_edge.end(), n) == new_edge.end())
                new_edge.push_back(n);
        }

        if (H.has_edge(new_edge))
            continue;
        
        H.add_edge(new_edge);

        return 1;
    }
    return 0;
}


// ============================================================================
// Loss functions
// ============================================================================

double degree_rmse(const Hypergraph &A, const Hypergraph &B)
{
    auto hist_A {A.degree_histogram()};
    auto hist_B {B.degree_histogram()};

    size_t max_size {std::max(hist_A.size(), hist_B.size())};
    hist_A.resize(max_size);
    hist_B.resize(max_size);

    double loss {};
    for (size_t i {}; i < max_size; ++i)
        loss += (hist_A[i] - hist_B[i]) * (hist_A[i] - hist_B[i]) / static_cast<double>(max_size);

    return std::sqrt(loss);
}

double degree_emd(const Hypergraph &A, const Hypergraph &B)
{
    auto hist_A {A.degree_histogram()};
    auto hist_B {B.degree_histogram()};

    double sum1 = 0.0;
    double sum2 = 0.0;

    for (double x : hist_A) { sum1 += x; }
    for (double x : hist_B) { sum2 += x; }

    const size_t n = std::max(hist_A.size(), hist_B.size());

    double cdf1 = 0.0;
    double cdf2 = 0.0;
    double distance = 0.0;

    for (size_t k = 0; k + 1 < n; ++k) {
        if (k < hist_A.size()) { cdf1 += hist_A[k] / sum1; }
        if (k < hist_B.size()) { cdf2 += hist_B[k] / sum2; }
        distance += std::abs(cdf1 - cdf2);
    }

    return distance;
}

double m3_difference(const Hypergraph& A, const std::vector<long long>& motif_target)
{
    auto motif_A {A.motif_count(3)};

    double loss {};
    for (size_t i {}; i < motif_A.size(); ++i) {
        double d {static_cast<double>(motif_A[i] - motif_target[i]) / (static_cast<double>(std::min(motif_target[i], motif_A[i])) + 1)};
        loss += d * d;
    }

    return std::sqrt(loss / static_cast<double>(motif_A.size()));
}

double m4_difference(const Hypergraph& A, const std::vector<long long>& motif_target, const std::vector<unsigned int>& motifs)
{
    auto motif_A {A.motif_count(4)};

    double loss {};
    for (auto m : motifs) {
        double d {static_cast<double>(motif_A[m] - motif_target[m]) / (static_cast<double>(std::min(motif_target[m], motif_A[m])) + 1)};
        loss += d * d;
    }

    return std::sqrt(loss / static_cast<double>(motifs.size()));
}


// ============================================================================
// Actual generator
// ============================================================================

double P_new_node(size_t N, size_t Nt, size_t M, size_t Mt)
{
    if (N >= Nt) { return 0.0; }

    double fn {static_cast<double>(N)/static_cast<double>(Nt)};
    double fm {static_cast<double>(M)/static_cast<double>(Mt)};

    double p {std::pow((1.0 - fn), 0.25) * (fm/fn)};

    return std::min(1.0, p);
}

size_t get_sampled_size(const std::vector<int>& size_hist)
{
    std::vector<int> sizes (size_hist.size());
    std::iota(sizes.begin(), sizes.end(), 0);
    return static_cast<size_t>(weighted_sample(sizes, size_hist));
}

std::vector<int> adapted_size_histogram(const std::vector<int>& target_size_hist, const std::vector<int>& gen_size_hist)
{
    std::vector<int> adapted_hist (target_size_hist.size());
    for (size_t i {}; i < target_size_hist.size(); ++i)
    {
        if (i < gen_size_hist.size())
            adapted_hist[i] = std::max(0, target_size_hist[i] - gen_size_hist[i]);
        else
            adapted_hist[i] = target_size_hist[i];
    }
    return adapted_hist;
}


Hypergraph generate_hypergraph(const std::vector<double>& weights, size_t Nt, size_t Mt, std::vector<int> size_hist, bool dynamic_weights, int max_attempts)
{
    unsigned int n_mech {5};
    static const std::vector<unsigned int> mechanism_ids {0, 1, 2, 3, 4};

    Hypergraph H {};

    ID initial_nodes {static_cast<ID>(std::ceil(0.05*static_cast<double>(Nt)))};
    for (ID n {}; n < initial_nodes; ++n) { H.add_node(n); }

    for (int i {}; i < static_cast<int>(std::ceil(0.05*static_cast<double>(Mt))); ++i) {
        random_edge(H, get_sampled_size(adapted_size_histogram(size_hist, H.size_histogram())), 0, max_attempts);
    }

    size_t i {};
    size_t N {H.node_count()};
    size_t M {H.edge_count()};

    std::vector<double> final_weights (n_mech);
    if (!dynamic_weights) { final_weights = weights; }
    
    while (M < Mt)
    {
        if (i > 10*Mt) { break; }

        ++i;

        if (dynamic_weights) {
            for (size_t i {}; i < n_mech; ++i) {
                double x {static_cast<double>(M)/static_cast<double>(Mt)};
                final_weights[i] = (1.0 - x)*weights[i] + x*weights[i+n_mech];
            }
        }

        auto g {weighted_sample(mechanism_ids, final_weights)};
        double p_new {P_new_node(N, Nt, M, Mt)};

        switch (g)
        {
        case 0:
            random_edge(H, get_sampled_size(adapted_size_histogram(size_hist, H.size_histogram())), p_new, max_attempts);
            break;
        case 1:
            PA(H, get_sampled_size(adapted_size_histogram(size_hist, H.size_histogram())), p_new, max_attempts);
            break;
        case 2:
            mutation(H, p_new, max_attempts);
            break;
        case 3:
            nesting(H, max_attempts);
            break;
        case 4:
            merging(H, max_attempts);
            break;
        }
        N = H.node_count();
        M = H.edge_count();
    }

    return H;
}


double test_weights(const Hypergraph& target, std::vector<double> weights, int loss_id, size_t n_hypergraphs, std::vector<unsigned int> motifs, bool two_phase, int max_attempts)
{
    auto N {target.node_count()};
    auto M {target.edge_count()};
    auto size_hist {target.size_histogram()};

    double loss {};

    std::vector<long long> motif_target {};
    if (loss_id == 3) { motif_target = target.motif_count(3); }

    else if (loss_id == 4) {
        motif_target = target.motif_count(4); 
        if (motifs.empty()) {
            motifs.resize(171);
            std::iota(motifs.begin(), motifs.end(), 0);
        }
    }
    
    for (size_t i {}; i < n_hypergraphs; ++i)
    {
        Hypergraph H = generate_hypergraph(weights, N, M, size_hist, two_phase, max_attempts);

        switch (loss_id)
        {
        case 1:
            loss += degree_emd(H, target);
            break;
        case 2:
            loss += degree_rmse(H, target);
            break;
        case 3:
            loss += m3_difference(H, motif_target);
            break;
        case 4:
            loss += m4_difference(H, motif_target, motifs);
            break;
        default:
            throw std::runtime_error("Unknown loss function.");
        }
    }

    return loss / static_cast<double>(n_hypergraphs);
}