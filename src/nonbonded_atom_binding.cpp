#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include "nonbonded_atom.cpp"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(nonbonded_atom, m) {
    nb::class_<Atom>(m, "Atom")
        .def(nb::init<>())
        .def("add_meta_properties", &Atom::add_meta_properties)
        .def("add_nonbonded_params", &Atom::add_nonbonded_params)
        .def("update_position", &Atom::update_position)
        .def("get_params", &Atom::get_params)
        .def("get_name", &Atom::get_name)
        .def("get_element", &Atom::get_element)
        .def("get_id", &Atom::get_id)
        .def("get_residue_id", &Atom::get_residue_id)
        .def("get_charge", &Atom::get_charge)
        .def("get_sigma", &Atom::get_sigma)
        .def("get_epsilon", &Atom::get_epsilon)
        .def("get_position", &Atom::get_position);
}