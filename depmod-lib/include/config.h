#pragma once

#include "stdafx.h"

class Configuration {
public:
  Configuration() = default;
  
  Configuration(double gammadot, double t_min, double t_max, size_t N, size_t K);

  double gammadot;
  double t_min;
  double t_max;
  double dt;
  size_t N;
  size_t K;
  size_t M;
};
