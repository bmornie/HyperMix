#include "hypergraph.hpp"

#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>

#include "unordered_dense.h"
#include "hypergraph_motif.hpp"


void Hypergraph::add_node(ID n)
{
    if (!m_nodes.contains(n))
        m_nodes[n] = {};
}

void Hypergraph::add_edge(std::vector<ID> nodes)
{
    std::sort(nodes.begin(), nodes.end());
    m_edges[m_next_id] = nodes;

    m_edge_ids.insert(nodes);

    for (ID n : nodes)
        m_nodes[n].push_back(m_next_id);

    ++m_next_id;
}


void Hypergraph::remove_edge(ID e)
{
    m_edge_ids.erase(m_edges[e]);

    for (ID n : m_edges[e]) {
        std::erase(m_nodes[n], e);
    }
    m_edges.erase(e);
}

void Hypergraph::add_node_to_edge(ID e, ID n)
{
    m_edge_ids.erase(m_edges[e]);

    auto pos = std::lower_bound(m_edges[e].begin(), m_edges[e].end(), n);
    m_edges[e].insert(pos, n);

    m_edge_ids.insert(m_edges[e]);
    
    m_nodes[n].push_back(e);
}

std::vector<ID> Hypergraph::nodes() const
{
    std::vector<ID> keys;
    keys.reserve(m_nodes.size());

    for (const auto &[key, _] : m_nodes)
        keys.push_back(key);

    return keys;
}

std::vector<ID> Hypergraph::edges() const
{
    std::vector<ID> keys;
    keys.reserve(m_edges.size());

    for (const auto &[key, _] : m_edges)
        keys.push_back(key);

    return keys;
}

std::vector<std::vector<ID>> Hypergraph::edge_list() const
{
    std::vector<std::vector<ID>> list {};
    list.reserve(m_edges.size());
    for (const auto &[key, value] : m_edges)
        list.push_back(value);
    return list;
}

int Hypergraph::degree(ID n, size_t edge_size) const
{
    if (!edge_size)
        return static_cast<int>(m_nodes.at(n).size());

    int deg{};
    for (auto e : m_nodes.at(n))
    {
        if (m_edges.at(e).size() == edge_size)
            ++deg;
    }
    return deg;
}

double Hypergraph::average_size() const
{
    double av_size{};
    int count{};
    for (const auto &[key, value] : m_edges)
    {
        av_size += static_cast<double>(value.size());
        ++count;
    }

    return av_size / count;
}

ankerl::unordered_dense::set<ID> Hypergraph::neighbors(ID n, size_t edge_size) const
{
    ankerl::unordered_dense::set<ID> nbs {};
    for (auto e : m_nodes.at(n))
    {
        if (!edge_size || edge_size == size(e))
        {
            for (auto nb : m_edges.at(e))
            {
                if (nb != n)
                    nbs.insert(nb);
            }
        }
    }
    return nbs;
}

std::vector<int> Hypergraph::degree_histogram() const
{
    if (m_nodes.empty()) { return {}; }

    std::vector<size_t> degrees (m_nodes.size());
    size_t i {};
    for (const auto &[n, edges] : m_nodes)
        degrees[i++] = edges.size();

    std::vector<int> histogram (*std::max_element(degrees.begin(), degrees.end()) + 1);
    for (auto d : degrees)
        histogram[d] += 1;

    return histogram;
}

std::vector<int> Hypergraph::size_histogram() const
{
    if (m_edges.empty()) { return {}; }

    std::vector<size_t> sizes (m_edges.size());
    size_t i {};
    for (const auto &[e, nodes] : m_edges)
        sizes[i++] = nodes.size();

    std::vector<int> histogram (*std::max_element(sizes.begin(), sizes.end()) + 1);
    for (auto s : sizes)
        histogram[s] += 1;

    return histogram;
}

std::vector<int> Hypergraph::intersection_histogram() const
{
    ankerl::unordered_dense::map<uint64_t, unsigned int> intersections {};

    for (const auto& [_, edge_ids] : m_nodes) {

        for (size_t i = 0; i < edge_ids.size(); ++i) {
            for (size_t j = i+1; j < edge_ids.size(); ++j) {

                auto e1 {edge_ids[i]};
                auto e2 {edge_ids[j]};
                if (e1 > e2) { std::swap(e1, e2); }

                uint64_t key {(uint64_t(e1) << 32) | e2 };
                ++ intersections[key];
            }
        }
    }

    std::vector<int> histogram {};
    for (const auto& [_, intersect] : intersections) {
        while (histogram.size() < static_cast<size_t>(intersect)) { histogram.push_back(0); }
        ++histogram[intersect-1];
    }
    return histogram;
}

std::vector<long long> Hypergraph::motif_count(int order) const 
{
    if (order == 3){
        Hypergraph_motif H3 {};
        for (const auto &[_, edge] : m_edges){
            if (edge.size() <= 3)
                H3.add_edge(edge);
        }
        return H3.motif_count_3();
    }

    if (order == 4){
        std::vector<long long> counts (171);

        Hypergraph_motif H4 {};
        for (const auto &[_, edge] : m_edges){
            if (edge.size() == 2)
                H4.add_edge(edge);
        }
        H4._motif_count_4_size_2(counts);

        for (const auto &[_, edge] : m_edges){
            if (edge.size() == 3)
                H4.add_edge(edge);
        }
        H4._motif_count_4_size_3(counts);

        for (const auto &[_, edge] : m_edges){
            if (edge.size() == 4)
                H4.add_edge(edge);
        }
        H4._motif_count_4_size_4(counts);

        return counts;
    }

    throw std::runtime_error("Computation of motif counts of order > 4 has not been implemented yet.");
}

Hypergraph read_hyperedge_list_from_stream(std::istream& input)
{
    Hypergraph H {};
    std::string line;

    while (std::getline(input, line)) {
        if (line.empty()) continue;

        std::replace(line.begin(), line.end(), ',', ' ');

        std::istringstream iss(line);
        std::vector<ID> edge;
        ID node;

        while (iss >> node) {
            if (std::find(edge.begin(), edge.end(), node) == edge.end()) {
                edge.push_back(node);
            }
        }

        if (!edge.empty()) {
            H.add_edge(edge);
        }
    }

    return H;
}

Hypergraph read_hyperedge_list(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    return read_hyperedge_list_from_stream(file);
}

Hypergraph read_hypergraph_string(const std::string& contents)
{
    std::istringstream input(contents);
    return read_hyperedge_list_from_stream(input);
}


std::vector<std::vector<size_t>> clique_expansion(const Hypergraph& H)
{
    const auto& hyperedges {H.edge_list()};
    std::unordered_map<ID, size_t> id_to_index {};

    for (const auto& edge : hyperedges) {
        for (ID id : edge) {
            id_to_index.emplace(id, id_to_index.size());
        }
    }

    std::vector<std::vector<size_t>> adj (id_to_index.size());

    for (const auto& edge : hyperedges) {
        for (size_t i = 0; i < edge.size(); ++i) {
            for (size_t j = i + 1; j < edge.size(); ++j) {
                size_t u = id_to_index.at(edge[i]);
                size_t v = id_to_index.at(edge[j]);

                adj[u].push_back(v);
                adj[v].push_back(u);
            }
        }
    }

    return adj;
}

double effective_diameter(const Hypergraph& H)
{
    std::vector<std::vector<size_t>> adj {clique_expansion(H)};
    const size_t n = adj.size();
    if (n < 2) return 0.0;

    std::vector<long long> hist(n, 0);
    std::vector<long long> dist(n);
    std::queue<size_t> q;

    long long total_reachable = 0;

    for (size_t s = 0; s < n; ++s) {
        std::fill(dist.begin(), dist.end(), -1);

        while (!q.empty()) q.pop();

        dist[s] = 0;
        q.push(s);

        while (!q.empty()) {
            size_t u = q.front();
            q.pop();

            for (size_t v : adj[u]) {
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    q.push(v);
                }
            }
        }

        // Count each reachable unordered pair exactly once.
        for (size_t v = s + 1; v < n; ++v) {
            if (dist[v] != -1) {
                ++hist[static_cast<size_t>(dist[v])];
                ++total_reachable;
            }
        }
    }

    if (total_reachable == 0) return 0.0;

    const double target = 0.9;
    long long cumulative = 0;
    double previous_fraction = 0.0;

    for (size_t d = 1; d < n; ++d) {
        cumulative += hist[d];

        const double fraction = static_cast<double>(cumulative) / static_cast<double>(total_reachable);

        if (fraction >= target) {
            const double increment = fraction - previous_fraction;

            if (increment == 0.0) {
                return static_cast<double>(d);
            }

            return static_cast<double>(d - 1)
                + (target - previous_fraction) / increment;
        }

        previous_fraction = fraction;
    }

    throw std::runtime_error(
        "Failed to compute effective diameter"
    );
}
