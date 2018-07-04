/*################################################################################
  ##
  ##   Copyright (C) 2017-2018 Odran Bonnet, 
  ##                           Alfred Galichon, 
  ##                           Keith O'Hara, and
  ##                           Matt Shum
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

#include "Rgeogram.hpp"
using namespace Rcpp;

int OTM2D(const double* points_inp, const int n_vert, const int vert_dim, double* weights_ptr);
int OTM3D(const double* points_inp, const int n_vert, const int vert_dim, double* weights_ptr);

SEXP OTM2D_R(SEXP eps_mat_R)
{
    try {
        int nrow = Rf_nrows(eps_mat_R);
        int ncol = Rf_ncols(eps_mat_R);

        if (ncol > nrow)
        {
            ::Rf_error( "Rgeogram: nrow < ncol!" );
            return R_NilValue;
        }

        //

        arma::vec weights_out = arma::zeros(nrow,1);

        Rgeogram::clocktime_t start_time = Rgeogram::tic();
        OTM2D(REAL(eps_mat_R), nrow, ncol, weights_out.memptr());
        Rgeogram::comptime_t algo_runtime = Rgeogram::tic() - start_time;

        double runtime_out = algo_runtime.count();

        //

        return Rcpp::List::create(Rcpp::Named("weights") = weights_out,
                                  Rcpp::Named("elapsed_time") = runtime_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "Rgeogram: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP OTM3D_R(SEXP eps_mat_R)
{
    try {
        int nrow = Rf_nrows(eps_mat_R);
        int ncol = Rf_ncols(eps_mat_R);

        if (ncol > nrow)
        {
            ::Rf_error( "Rgeogram: nrow < ncol!" );
            return R_NilValue;
        }

        //

        arma::vec weights_out = arma::zeros(nrow,1);

        Rgeogram::clocktime_t start_time = Rgeogram::tic();
        OTM3D(REAL(eps_mat_R), nrow, ncol, weights_out.memptr());
        Rgeogram::comptime_t algo_runtime = Rgeogram::tic() - start_time;

        double runtime_out = algo_runtime.count();

        //

        return Rcpp::List::create(Rcpp::Named("weights") = weights_out,
                                  Rcpp::Named("elapsed_time") = runtime_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "Rgeogram: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
