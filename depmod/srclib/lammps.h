#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

#include "config.h"
#include "deformation.h"

namespace py = pybind11;


void lammps_generate_box_evolution_data_brute(
  Configuration*, 
  Deformation*,
  py::array_t<double>,
  py::array_t<double>,
  py::array_t<double>,
  py::array_t<double>
);


void lmp_box_gen_brute_npts(
  Configuration*, 
  Deformation*,
  py::array_t<double>,
  py::array_t<double>,
  py::array_t<double>,
  py::array_t<double>
);


void lmp_box_gen_brute_kpts(
  Configuration*, 
  Deformation*,
  py::array_t<double>,
  py::array_t<double>,
  py::array_t<double>,
  py::array_t<double>
);


// inline void fill_array(const size_t index, const double dt, const double tfac, py::array_t<double>& out, const py::array_t<double>& G, const py::array_t<double>& F) {
//   auto o = out.mutable_unchecked<2>();
//   auto g = G.unchecked<2>();
//   auto f = F.unchecked<2>();

//   o(index, 0) = index;
//   o(index, 1) = index * dt * tfac;

//   o(index, 2) = g(index, 0); // lx
//   o(index, 3) = g(index, 4); // ly
//   o(index, 4) = g(index, 8); // lz
//   o(index, 5) = g(index, 5); // yz
//   o(index, 6) = g(index, 2); // xz
//   o(index, 7) = g(index, 1); // xy
//   
//   o(index,  8) = f(index, 0); // F_11
//   o(index,  9) = f(index, 1); // F_12
//   o(index, 10) = f(index, 2); // F_13
//   o(index, 11) = f(index, 3); // F_21
//   o(index, 12) = f(index, 4); // F_22
//   o(index, 13) = f(index, 5); // F_23
//   o(index, 14) = f(index, 6); // F_31
//   o(index, 15) = f(index, 7); // F_32
//   o(index, 16) = f(index, 8); // F_33
// };




// inline void lammps_generate_box_data(
//   size_t ntfs,
//   double dt,
//   double gdot,
//   double tfac,
//   py::array_t<double> Ht,
//   py::array_t<double> Gt,
//   py::array_t<double> Ft,
//   py::array_t<double> out,
//   Deformation* deform
// ) {

//   auto r_Ht = Ht.mutable_unchecked<2>();
//   auto r_Gt = Gt.mutable_unchecked<2>();
//   auto r_Ft = Ft.mutable_unchecked<2>();

//   // Const map for H0
//   Eigen::Map<const Matrix3d> H0(r_Ht.data(0, 0), 3, 3);
//   // Map for H(t)
//   Eigen::Map<Matrix3d> Hi(r_Ht.mutable_data(0, 0), 3, 3);
//   // Map for G(t)
//   Eigen::Map<Matrix3d> Gi(r_Gt.mutable_data(0, 0), 3, 3);
//   // Map for F(t)
//   Eigen::Map<Matrix3d> Fi(r_Ft.mutable_data(0, 0), 3, 3);
//   // Map for F(t-dt)
//   Eigen::Map<const Matrix3d> Fp(r_Ft.data(0, 0), 3, 3);
//   // Map for Row of out

//   // Handle special case i = 0 because we do start the loop at 1  
//   fill_array(0, dt, tfac, out, Gt, Ft);

//   // Loop to compute F increment
//   for (size_t i = 1; i < ntfs+1; i++) {
//     
//     // In principle data are swapped with no memory allocation
//     new (&Hi) Eigen::Map<Matrix3d>(r_Ht.mutable_data(i, 0), 3, 3);
//     new (&Gi) Eigen::Map<Matrix3d>(r_Gt.mutable_data(i, 0), 3, 3);
//     new (&Fi) Eigen::Map<Matrix3d>(r_Ft.mutable_data(i, 0), 3, 3);
//     
//     // F(t-dt)
//     new (&Fp) Eigen::Map<const Matrix3d>(r_Ft.data(i-1, 0), 3, 3);

//     // Compute F(t) = F(t-1) + dF(t)
//     Fi = Fp.array() + deform->get_Fdot_inc(dt, gdot, Fp, deform->S()).array();

//     // Compute the lattice H(t) = F(t)*H(0)
//     Hi = Fi * H0;

//     // Rotate back into lammps conventions (i.e. upper triangular matrix)
//     align_to_lammps_convention(Gi, Hi);

//     fill_array(i, dt, tfac, out, Gt, Ft);
//   };
// }
