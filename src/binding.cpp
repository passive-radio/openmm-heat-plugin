#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include "pairs.cpp"
// #include "nonbonded_atom.cpp"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(openmm_heat_plugin, m) {
    nb::class_<Pair>(m, "Pair")
        .def(nb::init<int, int>())
        .def(nb::init<int, int>())
        .def(nb::init<int, int, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>())
        .def("get_fij", &Pair::get_fij);

    nb::class_<Pairs>(m, "Pairs")
        .def(nb::init<int>())
        .def("add_pair", &Pairs::add_pair)
        .def("add_pairs", &Pairs::add_pairs)
        .def("reset", &Pairs::reset)
        .def("num_pairs", &Pairs::num_pairs);
}
