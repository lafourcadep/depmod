#pragma once

#include <Eigen/Dense>

template<typename T>
using MatrixABC = Eigen::MatrixBase<T>;


template<typename T, int R, int C>
using Matrix = Eigen::Matrix<T, R, C, Eigen::RowMajor>;


template<typename T>
using Ref = Eigen::Ref<T>;


template<typename T>
using Map = Eigen::Map<T>;


typedef Matrix<double, 3, 3> Matrix3d;
typedef Eigen::Vector3d Vector3d;
