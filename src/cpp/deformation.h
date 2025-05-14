#pragma once

#include "dtype.h"


class Deformation {
public:
  virtual ~Deformation() = default;
  // Precompute K = dt * gdot * S
  virtual Matrix3d compute_K(double, double) = 0;

  virtual const Matrix3d get_S(void) { return m_S; };

protected:
  Deformation(const Matrix3d& S);
  Matrix3d m_S;
};


class TractionCompression: public Deformation {
public:
  TractionCompression(const Matrix3d& S, bool cflag);
  TractionCompression(const Matrix3d& S);
  Matrix3d compute_K(double, double);

  inline bool cflag() { return m_cflag; }
  inline bool isoV() { return false; }

private:
  bool m_cflag;
};


class TractionCompressionIsoV: public Deformation {
public:
  TractionCompressionIsoV(const Matrix3d& S, bool cflag);
  TractionCompressionIsoV(const Matrix3d& S);
  Matrix3d compute_K(double, double);

  inline bool cflag() { return m_cflag; }
  inline bool isoV() { return true; }

private:
  bool m_cflag;
};


class PureShear: public Deformation {
public:
  PureShear(const Matrix3d& S);
  Matrix3d compute_K(double, double);
};


// class Deformation {
// public:
//   
//   virtual ~Deformation() = default;
//   virtual const Matrix3d get_Fdot_inc(double, double, const Matrix3d&, const Matrix3d&) = 0;

//   // const Matrix3d& S() { return m_S; };

//   const Matrix3d S;

// protected:
//   Deformation(const Matrix3d& S): m_S(S) {};
//   Matrix3d m_S;
// };

// class PureShear: public Deformation {
// public:
//   PureShear(const Matrix3d& S): Deformation(S) {};
//   const Matrix3d get_Fdot_inc(double dt, double gdot, const Matrix3d& F, const Matrix3d& S) override {
//     return (dt * gdot * (S * F).array()).matrix();
//   };
// };



// class Deformation {
// public:
//   virtual ~Deformation() {};
//   virtual int get_delta_Fdot(double, double, double) = 0;

// protected:
//   Eigen::Map<Mat33> m_S;
// };

// class PureShear: public Deformation {
// public:

//   PureShear(py::array_t<double> S) {
//     auto S_ = S.unchecked<2>();
//     m_S = Eigen::Map<const Mat33>(S_.data(0, 0), 3, 3);
//   };

//   int get_delta_Fdot(double, double, double);
// };

// class Traction: public Deformation {
// public:
//   int get_delta_Fdot(double, double, double);
// };

// class Compression: public Deformation {
// public:
//   int get_delta_Fdot(double, double, double);
// };

// class PyDeformation: public Deformation {
// public:
//   // Inherit the constructors
//   using Deformation::Deformation;

//   // Trampolines (need one for each virtual function)
//   int get_dFdot() {
//     PYBIND
//   }

// }
