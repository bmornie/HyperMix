#pragma once

#include <vector>
#include "hypergraph.hpp"



// ============================================================================
// Hyperedge formation mechanisms
// ============================================================================

int random_edge(Hypergraph& H, size_t edge_size, double p_new, int max_attempts);

int PA(Hypergraph& H, size_t edge_size, double p_new, int max_attempts);

int mutation(Hypergraph &H, double p_new, int max_attempts);

int nesting(Hypergraph& H, int max_attempts);

int merging(Hypergraph& H, int max_attempts);


// ============================================================================
// Loss functions
// ============================================================================

double degree_rmse(const Hypergraph &A, const Hypergraph &B);

double degree_emd(const Hypergraph &A, const Hypergraph &B);

double m3_difference(const Hypergraph& A, const std::vector<long long>& motif_target);

double m4_difference(const Hypergraph& A, const std::vector<long long>& motif_target, const std::vector<unsigned int>& motifs);


// ============================================================================
// Actual generator
// ============================================================================

double P_new_node(size_t N, size_t Nt, size_t M, size_t Mt);

size_t get_sampled_size(const std::vector<int>& size_hist);

std::vector<int> adapted_size_histogram(const std::vector<int>& target_size_hist, const std::vector<int>& gen_size_hist);

Hypergraph generate_hypergraph(const std::vector<double>& weights, size_t Nt, size_t Mt, std::vector<int> size_hist, bool dynamic_weights=true, int max_attempts=100);

double test_weights(const Hypergraph& target, std::vector<double> weights, int loss_id, size_t n_hypergraphs, std::vector<unsigned int> motifs=std::vector<unsigned int>{}, bool dynamic_weights=true, int max_attempts=100);