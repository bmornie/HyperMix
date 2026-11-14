#pragma once

#include <vector>
#include <string>

#include "unordered_dense.h"
#include "hypergraph_motif.hpp"


struct VectorHash
{
	size_t operator()(const std::vector<ID> &v) const
	{
		size_t hash = v.size();
		for (auto &i : v)
		{
			hash ^= std::hash<ID>{}(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		}
		return hash;
	}
};


class Hypergraph
{
private:
	ankerl::unordered_dense::map<ID, std::vector<ID>> m_nodes {};
	ankerl::unordered_dense::map<ID, std::vector<ID>> m_edges {};

	ankerl::unordered_dense::set<std::vector<ID>, VectorHash> m_edge_ids {};
	ID m_next_id {};


public:
	void add_node(ID n);

	void add_edge(std::vector<ID> nodes);

	std::vector<ID> get_edge(ID e) const { return m_edges.at(e); }

	void remove_edge(ID e);

	void add_node_to_edge(ID e, ID n);

	std::vector<ID> nodes() const;

	std::vector<ID> edges() const;

	size_t node_count() const { return m_nodes.size(); }

	size_t edge_count() const { return m_edges.size(); }

	bool has_edge(std::vector<ID> nodes) const
	{
		std::sort(nodes.begin(), nodes.end());
		return m_edge_ids.contains(nodes);
	}

	std::vector<std::vector<ID>> edge_list() const;

	int degree(ID n, size_t edge_size = 0) const;

	size_t size(ID e) const { return m_edges.at(e).size(); }

	double average_size() const;

	ankerl::unordered_dense::set<ID> neighbors(ID n, size_t edge_size = 0) const;

	std::vector<int> degree_histogram() const;

	std::vector<int> size_histogram() const;

	std::vector<int> intersection_histogram() const;

	std::vector<long long> motif_count(int order = 3) const;
};


// Core parser: reads from any input stream
Hypergraph read_hyperedge_list_from_stream(std::istream& input);

// Reads from a filepath
Hypergraph read_hyperedge_list(const std::string& filename);

// Reads from an in-memory string
Hypergraph read_hypergraph_string(const std::string& contents);


// ============================================================================
// Effective diameter
// ============================================================================

// Convert hyperedges to the adjacency list of the clique expansion.
std::vector<std::vector<size_t>> clique_expansion(const Hypergraph& H);

// 90th-percentile effective diameter with interpolation.
double effective_diameter(const Hypergraph& H);