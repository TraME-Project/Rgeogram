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

#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>

//

#include <exploragram/optimal_transport/optimal_transport.h>
#include <exploragram/optimal_transport/linear_least_squares.h>
#include <exploragram/optimal_transport/sampling.h>
#include <exploragram/optimal_transport/optimal_transport_3d.h>
#include <exploragram/optimal_transport/optimal_transport_2d.h>
#include <exploragram/optimal_transport/optimal_transport_on_surface.h>

#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_sampling.h>
#include <geogram/mesh/mesh_AABB.h>

#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/generic_RVD_vertex.h>
#include <geogram/voronoi/RVD_callback.h>
#include <geogram/voronoi/generic_RVD_cell.h>
#include <geogram/delaunay/delaunay_3d.h>

#include <geogram/points/nn_search.h>
#include <geogram/numerics/optimizer.h>
#include <geogram/numerics/lbfgs_optimizers.h>

#include <geogram/basic/stopwatch.h>
#include <geogram/basic/file_system.h>
#include <geogram/basic/process.h>
#include <geogram/basic/progress.h>
#include <geogram/basic/permutation.h>
#include <geogram/basic/command_line.h>

#include <geogram/NL/nl.h>

//

using uint_t = unsigned int;

namespace GEO
{

inline
void
set_mesh_point(Mesh &M, index_t v, const double *coords, index_t dim)
{
    geo_debug_assert(M.vertices.dimension() >= dim);
    double *p = M.vertices.point_ptr(v);
    for (index_t c = 0; c < dim; ++c)
    {
        p[c] = coords[c];
    }
}

}

using namespace GEO;

//
// 2D

int
OTM2D(const double* points_inp, const int n_vert, const int vert_dim, double* weights_ptr)
{
    
    // Initialize the Geogram library.
    GEO::initialize();

    // Import standard command line arguments.
    GEO::CmdLine::import_arg_group("standard");

    //
    GEO::CmdLine::set_arg("algo:predicates", "exact");

    //
    // Declare a mesh for sample points

    GEO::Mesh M_points;

    M_points.vertices.create_vertices(n_vert);

    double xyz[3];

    for (uint_t i=0U; i < n_vert; i++)
    {
        for (uint_t k=0U; k < vert_dim; k++)
        {
            xyz[k] = points_inp[i + k*n_vert]; // points_inp are in column-stacked order
        }

        GEO::set_mesh_point(M_points,i,xyz,vert_dim);
    }

    //
    // create a flat mesh for the uniform density

    GEO::Mesh M_unif;

    M_unif.vertices.create_vertices(4);

    xyz[0] = 1.0; xyz[1] = 1.0; xyz[2] = 0.0;
    GEO::set_mesh_point(M_unif,0,xyz,3);

    xyz[0] = 1.0; xyz[1] = 0.0; xyz[2] = 0.0;
    GEO::set_mesh_point(M_unif,1,xyz,3);

    xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
    GEO::set_mesh_point(M_unif,2,xyz,3);

    xyz[0] = 0.0; xyz[1] = 1.0; xyz[2] = 0.0;
    GEO::set_mesh_point(M_unif,3,xyz,3);

    // create facets

    vector<index_t> facet_vertices;
    facet_vertices.resize(0);
    facet_vertices.push_back(2); facet_vertices.push_back(1); facet_vertices.push_back(0);

    index_t f = M_unif.facets.create_polygon(facet_vertices.size());
    for (index_t lv = 0; lv < facet_vertices.size(); ++lv)
    {
        M_unif.facets.set_vertex(f, lv, facet_vertices[lv]);
    }

    facet_vertices.resize(0);
    facet_vertices.push_back(2); facet_vertices.push_back(0); facet_vertices.push_back(3);

    f = M_unif.facets.create_polygon(facet_vertices.size());
    for (index_t lv = 0; lv < facet_vertices.size(); ++lv)
    {
        M_unif.facets.set_vertex(f, lv, facet_vertices[lv]);
    }

    //
    // build the optimal transport map

    GEO::OptimalTransportMap2d OTM(&M_unif); 

    OTM.set_points(M_points.vertices.nb(),
                   M_points.vertices.point_ptr(0),
                   M_points.vertices.dimension());

    for (int i=0; i < M_points.vertices.nb(); i++)
    {
        OTM.set_nu(i, 1.0/double(M_points.vertices.nb()));
    }

    OTM.set_verbose(false);
    OTM.set_regularization(0.0);
    OTM.set_epsilon(1e-6);
    OTM.set_Newton(true);
    OTM.optimize(1000);

    // GEO::Mesh laguerre; // if we wanted the Laguerre tessellation

    // laguerre.clear();
    // OTM.get_RVD(laguerre);
    //

    // normalize the potentials by weight(0)

    // weights.set_size(OTM.nb_points());

    double weight0 = OTM.weight(0);
    for (int jj=0; jj < OTM.nb_points(); jj++)
    {
        weights_ptr[jj] = OTM.weight(jj) - weight0;
    }

    return 0;
}

//
// 3D

int
OTM3D(const double* points_inp, const int n_vert, const int vert_dim, double* weights_ptr)
{
    
    // Initialize the Geogram library.
    GEO::initialize();

    // Import standard command line arguments.
    GEO::CmdLine::import_arg_group("standard");

    //
    GEO::CmdLine::set_arg("algo:predicates", "exact");

    //
    // Declare a mesh for sample points

    GEO::Mesh M_points;

    M_points.vertices.create_vertices(n_vert);

    double xyz[3];

    for (uint_t i=0U; i < n_vert; i++)
    {
        for (uint_t k=0U; k < vert_dim; k++)
        {
            xyz[k] = points_inp[i + k*n_vert]; // points_inp are in column-stacked order
        }

        GEO::set_mesh_point(M_points,i,xyz,vert_dim);
    }

    //
    // create a flat mesh for the uniform density

    GEO::Mesh M_unif;

    M_unif.vertices.create_vertices(9);

    xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0;
    GEO::set_mesh_point(M_unif,0,xyz,3);

    xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 1.0;
    GEO::set_mesh_point(M_unif,1,xyz,3);

    xyz[0] = 0.0; xyz[1] = 1.0; xyz[2] = 0.0;
    GEO::set_mesh_point(M_unif,2,xyz,3);

    xyz[0] = 0.0; xyz[1] = 1.0; xyz[2] = 1.0;
    GEO::set_mesh_point(M_unif,3,xyz,3);

    xyz[0] = 1.0; xyz[1] = 0.0; xyz[2] = 0.0;
    GEO::set_mesh_point(M_unif,4,xyz,3);

    xyz[0] = 1.0; xyz[1] = 0.0; xyz[2] = 1.0;
    GEO::set_mesh_point(M_unif,5,xyz,3);

    xyz[0] = 1.0; xyz[1] = 1.0; xyz[2] = 0.0;
    GEO::set_mesh_point(M_unif,6,xyz,3);

    xyz[0] = 1.0; xyz[1] = 1.0; xyz[2] = 1.0;
    GEO::set_mesh_point(M_unif,7,xyz,3);

    xyz[0] = 0.494941; xyz[1] = 0.652018; xyz[2] = 0.319279;
    GEO::set_mesh_point(M_unif,8,xyz,3);

    // create cells

    M_unif.cells.create_tets(12);

    // for (index_t t = 0; t < nb_cells; ++t)
    // {
    //     for (index_t i = 0; i < 4; ++i)
    //     {
    //         index_t v = in.field_as_uint(i + 1);
    //         M_unif.cells.set_vertex(t, i, v);
    //     }
    // }

    M_unif.cells.set_vertex(0, 0, 0); M_unif.cells.set_vertex(0, 1, 2); 
    M_unif.cells.set_vertex(0, 2, 3); M_unif.cells.set_vertex(0, 3, 8);

    M_unif.cells.set_vertex(1, 0, 6); M_unif.cells.set_vertex(1, 1, 3); 
    M_unif.cells.set_vertex(1, 2, 8); M_unif.cells.set_vertex(1, 3, 7);

    M_unif.cells.set_vertex(2, 0, 1); M_unif.cells.set_vertex(2, 1, 8); 
    M_unif.cells.set_vertex(2, 2, 3); M_unif.cells.set_vertex(2, 3, 7);

    M_unif.cells.set_vertex(3, 0, 1); M_unif.cells.set_vertex(3, 1, 8); 
    M_unif.cells.set_vertex(3, 2, 5); M_unif.cells.set_vertex(3, 3, 0);

    M_unif.cells.set_vertex(4, 0, 8); M_unif.cells.set_vertex(4, 1, 1); 
    M_unif.cells.set_vertex(4, 2, 5); M_unif.cells.set_vertex(4, 3, 7);

    M_unif.cells.set_vertex(5, 0, 1); M_unif.cells.set_vertex(5, 1, 8); 
    M_unif.cells.set_vertex(5, 2, 0); M_unif.cells.set_vertex(5, 3, 3);

    M_unif.cells.set_vertex(6, 0, 5); M_unif.cells.set_vertex(6, 1, 8); 
    M_unif.cells.set_vertex(6, 2, 4); M_unif.cells.set_vertex(6, 3, 0);

    M_unif.cells.set_vertex(7, 0, 5); M_unif.cells.set_vertex(7, 1, 6); 
    M_unif.cells.set_vertex(7, 2, 8); M_unif.cells.set_vertex(7, 3, 7);

    M_unif.cells.set_vertex(8, 0, 5); M_unif.cells.set_vertex(8, 1, 6); 
    M_unif.cells.set_vertex(8, 2, 4); M_unif.cells.set_vertex(8, 3, 8);

    M_unif.cells.set_vertex(9, 0, 4); M_unif.cells.set_vertex(9, 1, 2); 
    M_unif.cells.set_vertex(9, 2, 0); M_unif.cells.set_vertex(9, 3, 8);

    M_unif.cells.set_vertex(10, 0, 6); M_unif.cells.set_vertex(10, 1, 3); 
    M_unif.cells.set_vertex(10, 2, 2); M_unif.cells.set_vertex(10, 3, 8);

    M_unif.cells.set_vertex(11, 0, 2); M_unif.cells.set_vertex(11, 1, 4); 
    M_unif.cells.set_vertex(11, 2, 6); M_unif.cells.set_vertex(11, 3, 8);

    //
    // build the optimal transport map

    M_unif.vertices.set_dimension(4);  // update the cube dimensions to dim + 1

    GEO::OptimalTransportMap3d OTM(&M_unif);

    OTM.set_points(M_points.vertices.nb(),
                   M_points.vertices.point_ptr(0),
                   M_points.vertices.dimension());

    for (int i=0; i < M_points.vertices.nb(); i++)
    {
        OTM.set_nu(i, 1.0/double(M_points.vertices.nb()));
    }

    OTM.set_verbose(false);
    OTM.set_regularization(0.0);
    OTM.set_epsilon(1e-6);
    OTM.set_Newton(true);
    OTM.optimize(1000);

    // GEO::Mesh laguerre; // if we wanted the Laguerre tessellation

    // laguerre.clear();
    // OTM.get_RVD(laguerre);
    //

    // normalize the potentials by weight(0)

    // weights.set_size(OTM.nb_points());

    double weight0 = OTM.weight(0);
    for (int jj=0; jj < OTM.nb_points(); jj++)
    {
        weights_ptr[jj] = OTM.weight(jj) - weight0;
    }

    return 0;
}
