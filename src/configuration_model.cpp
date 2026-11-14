#include "configuration_model.hpp"

#include <vector>
#include <algorithm>
#include <unordered_set>
#include <cmath>

#include "random.hpp"
#include "hypergraph.hpp"



void pairwise_shuffle(std::vector<ID>& edge1, std::vector<ID>& edge2)
{
    std::vector<ID> new1;
    std::vector<ID> new2;
    new1.reserve(edge1.size());
    new2.reserve(edge2.size());

    std::vector<ID> unique_nodes {};
    for (auto n : edge1){
        if (std::find(edge2.begin(), edge2.end(), n) == edge2.end()){
            unique_nodes.push_back(n);
        }
        else {
            new1.push_back(n);
            new2.push_back(n);
        }
    }
    for (auto n : edge2){
        if (std::find(edge1.begin(), edge1.end(), n) == edge1.end())
            unique_nodes.push_back(n);
    }

    while (unique_nodes.size() > 0){
        ID n {pop_random(unique_nodes)};
        if (new1.size() < edge1.size())
            new1.push_back(n);
        else
            new2.push_back(n);
    }

    std::sort(new1.begin(), new1.end());
    std::sort(new2.begin(), new2.end());

    edge1 = std::move(new1);
    edge2 = std::move(new2);
}

Hypergraph sample_configuration_model(const Hypergraph& H_orig, int steps)
{
    std::vector<std::vector<ID>> edge_list {H_orig.edge_list()};
    std::unordered_set<std::vector<ID>, VectorHash> edge_set {edge_list.begin(), edge_list.end()};
    auto m {edge_list.size()-1};

    int t {};
    while (t < steps){
        auto i {rand_int<std::size_t>(0, m)};
        auto j {rand_int<std::size_t>(0, m)};
        while (i == j)
            j = rand_int<std::size_t>(0, m);

        auto new_edge1 {edge_list[i]};
        auto new_edge2 {edge_list[j]};
        
        pairwise_shuffle(new_edge1, new_edge2);

        if (!edge_set.contains(new_edge1) && !edge_set.contains(new_edge2)){
            edge_set.erase(edge_list[i]);
            edge_set.erase(edge_list[j]);
            edge_set.insert(new_edge1);
            edge_set.insert(new_edge2);
            edge_list[i] = std::move(new_edge1);
            edge_list[j] = std::move(new_edge2);
            ++t;
        }
    }

    Hypergraph H {};
    for (const auto& edge : edge_list)
        H.add_edge(edge);
    
    return H;
}


std::vector<double> motif_zscore(const Hypergraph& H, int motif_size, int samples, int swaps)
{
    if (swaps == -1)
        swaps = 10*static_cast<int>(H.edge_count());
    
    auto motif_counts {H.motif_count(motif_size)};
    std::size_t m {motif_counts.size()};
    
    std::vector<double> sum (m);
    std::vector<double> sum_sq (m);

    for (int s {}; s < samples; ++s){
        Hypergraph Hs {sample_configuration_model(H, swaps)};
        auto counts {Hs.motif_count(motif_size)};

        for (std::size_t i {}; i < m; ++i){
            sum[i] += static_cast<double>(counts[i]);
            sum_sq[i] += static_cast<double>(counts[i]*counts[i]);
        }
    }

    std::vector<double> zscores (m);

    for (std::size_t i {}; i < m; ++i){
        double mean {sum[i]/static_cast<double>(samples)};
        double std {std::sqrt(sum_sq[i]/static_cast<double>(samples) - mean*mean)};
        zscores[i] = (static_cast<double>(motif_counts[i]) - mean) / std;
    }

    return zscores; 
}

std::vector<double> motif_abundance(const Hypergraph& H, int motif_size, int samples, int swaps, double epsilon)
{
    if (swaps == -1)
        swaps = 10*static_cast<int>(H.edge_count());
    
    auto motif_counts {H.motif_count(motif_size)};
    std::size_t m {motif_counts.size()};
    
    std::vector<long long> sum (m);

    for (int s {}; s < samples; ++s){
        Hypergraph Hs {sample_configuration_model(H, swaps)};
        auto counts {Hs.motif_count(motif_size)};

        for (std::size_t i {}; i < m; ++i){
            sum[i] += counts[i];
        }
    }

    std::vector<double> abundance (m);

    for (std::size_t i {}; i < m; ++i){
        double mean {static_cast<double>(sum[i])/static_cast<double>(samples)};
        abundance[i] = (static_cast<double>(motif_counts[i]) - mean) / (static_cast<double>(motif_counts[i]) + mean + epsilon);
    }

    return abundance; 
}

std::vector<std::vector<double>> motif_statistics(const Hypergraph& H, int motif_size, int samples, int swaps, double epsilon)
{
    if (swaps == -1)
        swaps = 10*static_cast<int>(H.edge_count());
    
    auto motif_counts {H.motif_count(motif_size)};
    std::size_t m {motif_counts.size()};
    
    std::vector<long long> sum (m);
    std::vector<double> sum_sq (m);

    for (int s {}; s < samples; ++s){
        Hypergraph Hs {sample_configuration_model(H, swaps)};
        auto counts {Hs.motif_count(motif_size)};

        for (std::size_t i {}; i < m; ++i){
            sum[i] += counts[i];
            sum_sq[i] += static_cast<double>(counts[i]*counts[i]);
        }
    }

    std::vector<double> zscores (m);
    std::vector<double> abundance (m);

    for (std::size_t i {}; i < m; ++i){
        double mean {sum[i]/static_cast<double>(samples)};
        double std {std::sqrt(sum_sq[i]/static_cast<double>(samples) - mean*mean)};
        zscores[i] = (static_cast<double>(motif_counts[i]) - mean) / std;
        abundance[i] = (static_cast<double>(motif_counts[i]) - mean) / (static_cast<double>(motif_counts[i]) + mean + epsilon);
    }

    return {zscores, abundance}; 
}