#pragma once

#include <cstdint>
#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "config.h"
#include "deformation.h"

namespace py = pybind11;

typedef std::function<void(Configuration*, Deformation*, py::array_t<double>, py::array_t<double>)> xsp_function_t;

const std::map<uint32_t, xsp_function_t>& exastamp_methods();

void _exastamp_box_gen_brute_n(Configuration*, Deformation*, py::array_t<double>, py::array_t<double>);
void _exastamp_box_gen_brute_k(Configuration*, Deformation*, py::array_t<double>, py::array_t<double>);
