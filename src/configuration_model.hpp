#pragma once

#include <vector>
#include "hypergraph.hpp"



void pairwise_shuffle(std::vector<ID>& edge1, std::vector<ID>& edge2);

Hypergraph sample_configuration_model(const Hypergraph& H_orig, int steps);


std::vector<double> motif_zscore(const Hypergraph& H, int motif_size, int samples=100, int swaps=-1);

std::vector<double> motif_abundance(const Hypergraph& H, int motif_size, int samples=100, int swaps=-1, double epsilon=4.0);

// Faster way to compute both zscores and abundance
std::vector<std::vector<double>> motif_statistics(const Hypergraph& H, int motif_size, int samples=100, int swaps=-1, double epsilon=4.0);