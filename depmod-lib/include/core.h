#pragma once

#include <cstdint>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "config.h"
#include "deformation.h"

namespace py = pybind11;

enum INTEGRATION_METHOD: uint32_t {
  BRUTE  = 1 << 0,
  CENTER = 1 << 1,
  BRUZY  = 1 << 2,
  KPTS   = 1 << 7
};

void lammps_generate_box_evolution_data(INTEGRATION_METHOD, Configuration*, Deformation*, py::array_t<double>, py::array_t<double>, py::array_t<double>, py::array_t<double>);

void exastamp_generate_box_evolution_data(INTEGRATION_METHOD, Configuration*, Deformation*, py::array_t<double>, py::array_t<double>);

// LAMMPS
// using lmp_func_t = std::function<void(Configuration*, Deformation*, py::array_t<double>, py::array_t<double>, py::array_t<double>, py::array_t<double>)>;

// struct lmp_method_t {
//   uint8_t flag;
//   lmp_func_t func;
// };

// const std::vector<lmp_method_t> lammps_methods();

// const lmp_func_t& get_lammps_method(uint32_t method_flag);


// void _lammps_box_gen_brute_n(Configuration*, Deformation*, py::array_t<double>, py::array_t<double>, py::array_t<double>, py::array_t<double>);
// void _lammps_box_gen_brute_k(Configuration*, Deformation*, py::array_t<double>, py::array_t<double>, py::array_t<double>, py::array_t<double>);
// void _lammps_box_gen_center_n();
// void _lammps_box_gen_center_k();
// void _lammps_box_gen_bruzy_n();
// void _lammps_box_gen_bruzy_k();


// EXASTAMP
// struct xsp_method_t {
//   uint8_t flag;
//   std::function<void(Configuration*, Deformation*, py::array_t<double>, py::array_t<double>)> func;
// };

// const std::vector<xsp_method_t> xstamp_methods();

// xsp_method_t& get_xstamp_method(uint32_t method_flag);


// void _xstamp_box_gen_brute_n();
// void _xstamp_box_gen_brute_k();
// void _xstamp_box_gen_center_n();
// void _xstamp_box_gen_center_k();
// void _xstamp_box_gen_bruzy_n();
// void _xstamp_box_gen_bruzy_k();
