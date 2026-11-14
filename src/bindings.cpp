#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "hypergraph.hpp"
#include "generator.hpp"
#include "configuration_model.hpp"

namespace py = pybind11;

PYBIND11_MODULE(hypermix, m) {
    py::class_<Hypergraph>(m, "Hypergraph")
        .def(py::init<>())
        .def("add_node", &Hypergraph::add_node)
        .def("add_edge", &Hypergraph::add_edge)
        .def("add_node_to_edge", &Hypergraph::add_node_to_edge)
        .def("nodes", &Hypergraph::nodes)
        .def("edges", &Hypergraph::edges)
        .def("node_count", &Hypergraph::node_count)
        .def("edge_count", &Hypergraph::edge_count)
        .def("get_edge", &Hypergraph::get_edge)
        .def("has_edge", &Hypergraph::has_edge)
        .def("edge_list", &Hypergraph::edge_list)
        .def("degree", &Hypergraph::degree)
        .def("size", &Hypergraph::size)
        .def("average_size", &Hypergraph::average_size)
        .def("neighbors", &Hypergraph::neighbors)
        .def("degree_histogram", &Hypergraph::degree_histogram)
        .def("size_histogram", &Hypergraph::size_histogram)
        .def("intersection_histogram", &Hypergraph::intersection_histogram)
        .def("motif_count", &Hypergraph::motif_count, py::arg("order")=3);
    
    m.def("read_hyperedge_list", &read_hyperedge_list, py::arg("filename"));

    m.def("read_hypergraph_string", &read_hypergraph_string, py::arg("contents"));

    m.def("effective_diameter", &effective_diameter, py::arg("H"));



    m.def("generate_hypergraph", &generate_hypergraph, py::arg("weights"), py::arg("Nt"), py::arg("Mt"), py::arg("size_hist"), py::arg("dynamic_weights")=true, py::arg("max_attempts")=100);

    m.def("test_weights", &test_weights, py::arg("target"), py::arg("weights"), py::arg("loss_id"), py::arg("n_hypergraphs"), py::arg("motifs")=std::vector<unsigned int>{}, py::arg("dynamic_weights")=true, py::arg("max_attempts")=100, py::call_guard<py::gil_scoped_release>());



    m.def("sample_configuration_model", &sample_configuration_model, py::arg("H_orig"), py::arg("steps"));

    m.def("motif_zscore", &motif_zscore, py::arg("H"), py::arg("motif_size"), py::arg("samples")=100, py::arg("swaps")=-1);

    m.def("motif_abundance", &motif_abundance, py::arg("H"), py::arg("motif_size"), py::arg("samples")=100, py::arg("swaps")=-1, py::arg("epsilon")=4.0);

    m.def("motif_statistics", &motif_statistics, py::arg("H"), py::arg("motif_size"), py::arg("samples")=100, py::arg("swaps")=-1, py::arg("epsilon")=4.0);
}