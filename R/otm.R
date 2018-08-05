################################################################################
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
################################################################################

otm2D <- function(chi_mat, weights_in = NULL)
{
    if (is.null(weights_in))
    {
        weights_in <- rep(1.0/nrow(chi_mat),nrow(chi_mat))
    }

    if (sum(weights_in) != 1)
    {
        weights_in <- weights_in / sum(weights_in)
    }

    res <- .Call("OTM2D_R", chi_mat, weights_in, PACKAGE = "Rgeogram")

    return(res)
}

otm3D <- function(chi_mat, weights_in = NULL)
{
    if (is.null(weights_in))
    {
        weights_in <- rep(1.0/nrow(chi_mat),nrow(chi_mat))
    }

    if (sum(weights_in) != 1)
    {
        weights_in <- weights_in / sum(weights_in)
    }

    #

    res <- .Call("OTM3D_R", chi_mat, weights_in, PACKAGE = "Rgeogram")

    return(res)
}
