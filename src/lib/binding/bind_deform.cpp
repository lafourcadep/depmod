#include "pybind11/cast.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

#include "deformation.h"

/* Trampoline class allowing for extension of Deformation class directly from Python */

class PyDeformation: public Deformation {
public:
  // Inherit constructors
  // using Deformation::Deformation;
  PyDeformation(const Matrix3d& S): Deformation(S) {};
  
  // Trampoline function (need one for each virtual function)
  Matrix3d compute_K(double dt, double gammadot) override {
    PYBIND11_OVERRIDE_PURE(
      Matrix3d, // return type
      Deformation, // parent class
      compute_K, // Function name
      dt, // input arguments
      gammadot
    );
  };
};

void init_deformation(py::module_& m) {
  
  auto subm = m.def_submodule("deformation");

  subm.doc() = "Deformation sub-module.";

  py::class_<Deformation, PyDeformation>(subm, "_lib_Deformation")
    .def(py::init<const Matrix3d&>())
    .def("compute_K", &Deformation::compute_K)
    .def_property_readonly("S", &Deformation::get_S);

  py::class_<TractionCompression, Deformation>(subm, "_lib_TractionCompression")
    .def(py::init<const Matrix3d&, bool>())
    .def(py::init<const Matrix3d&>())
    .def("flags", [](TractionCompression& a) { return py::make_tuple(a.cflag(), a.isoV()); });

  py::class_<TractionCompressionIsoV, Deformation>(subm, "_lib_TractionCompressionIsoV")
    .def(py::init<const Matrix3d&, bool>())
    .def(py::init<const Matrix3d&>())
    .def("flags", [](TractionCompressionIsoV& a) { return py::make_tuple(a.cflag(), a.isoV()); });

  py::class_<PureShear, Deformation>(subm, "_lib_PureShear")
    .def(py::init<const Matrix3d&>());
};
