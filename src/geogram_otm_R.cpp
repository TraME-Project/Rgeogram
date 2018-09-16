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

#include "Rgeogram.hpp"
using namespace Rcpp;

SEXP OTM2D_R(SEXP chi_mat_R, SEXP weights_in_R)
{
    try {
        Rgeogram::uint_t nrow = Rf_nrows(chi_mat_R);
        Rgeogram::uint_t ncol = Rf_ncols(chi_mat_R);

        if (ncol > nrow)
        {
            ::Rf_error( "Rgeogram: nrow < ncol!" );
            return R_NilValue;
        }

        //

        // arma::vec weights_out = arma::zeros(nrow,1);
        std::vector<double> weights_out(nrow,0);

        Rgeogram::clocktime_t start_time = Rgeogram::tic();
        Rgeogram::OTM2D(REAL(chi_mat_R), nrow, ncol, REAL(weights_in_R), weights_out.data());
        Rgeogram::comptime_t algo_runtime = Rgeogram::tic() - start_time;

        double runtime_out = algo_runtime.count();

        //

        return Rcpp::List::create(Rcpp::Named("weights") = wrap(weights_out),
                                  Rcpp::Named("elapsed_time") = runtime_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "Rgeogram: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP OTM3D_R(SEXP chi_mat_R, SEXP weights_in_R)
{
    try {
        Rgeogram::uint_t nrow = Rf_nrows(chi_mat_R);
        Rgeogram::uint_t ncol = Rf_ncols(chi_mat_R);

        if (ncol > nrow)
        {
            ::Rf_error( "Rgeogram: nrow < ncol!" );
            return R_NilValue;
        }

        //

        // arma::vec weights_out = arma::zeros(nrow,1);
        std::vector<double> weights_out(nrow,0);

        Rgeogram::clocktime_t start_time = Rgeogram::tic();
        Rgeogram::OTM3D(REAL(chi_mat_R), nrow, ncol, REAL(weights_in_R), weights_out.data());
        Rgeogram::comptime_t algo_runtime = Rgeogram::tic() - start_time;

        double runtime_out = algo_runtime.count();

        //

        return Rcpp::List::create(Rcpp::Named("weights") = wrap(weights_out),
                                  Rcpp::Named("elapsed_time") = runtime_out);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "Rgeogram: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
