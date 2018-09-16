/*################################################################################
  ##
  ##   Copyright (C) 2018 Keith O'Hara
  ##
  ##   This file is part of the Rgeogram package.
  ##
  ##   Rgeogram is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   Rgeogram is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with Rgeogram. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

#ifndef _RGEOGRAM_HPP
#define _RGEOGRAM_HPP

#include <Rcpp.h>
#include <chrono>

namespace Rgeogram
{
    using uint_t = unsigned int;

    #include "tictoc.hpp"

    int OTM2D(const double* points_inp, const uint_t n_vert, const uint_t vert_dim, const double* weights_in_ptr, double* weights_out_ptr);
    int OTM3D(const double* points_inp, const uint_t n_vert, const uint_t vert_dim, const double* weights_in_ptr, double* weights_out_ptr);
}

RcppExport SEXP OTM2D_R(SEXP chi_mat_R, SEXP weights_in_R);
RcppExport SEXP OTM3D_R(SEXP chi_mat_R, SEXP weights_in_R);

#endif
