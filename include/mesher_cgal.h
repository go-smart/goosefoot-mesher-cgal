/**
 * mesher_cgal
 *
 * Copyright (C) 2013-  NUMA Engineering Ltd. (see AUTHORs file)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MESHER_CGAL_H
#define MESHER_CGAL_H

#include "cgalsettings.pb.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>

// EPECK kernel
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// surface mesh
#include <CGAL/Polyhedron_3.h>

// nef
#include <CGAL/Nef_polyhedron_3.h> 
#include <CGAL/Nef_3/SNC_indexed_items.h>

#include <CGAL/bounding_box.h>


// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
typedef CGAL::Polyhedron_3<Exact_Kernel> Exact_polyhedron;
typedef CGAL::Cartesian_converter<Exact_Kernel,K> EK_to_K;
typedef CGAL::Cartesian_converter<K,Exact_Kernel> K_to_EK;
typedef CGAL::Polyhedral_mesh_domain_3<Exact_polyhedron, Exact_Kernel> Exact_Mesh_domain;
typedef CGAL::Nef_polyhedron_3<Exact_Kernel,
			       CGAL::SNC_indexed_items,
			       bool> Nef_polyhedron; 

typedef K::Point_3 Point;
typedef Exact_Kernel::Point_3 Exact_Point;
typedef K::FT FT;
typedef FT (Function)(const Point&);

typedef Exact_polyhedron::HalfedgeDS Exact_HalfedgeDS;

namespace mesherCGAL {
    //int visualization_set_allocation_order(int ix, int order, double cl, double needle_dist, double x, double y, double z);
    //int visualization_create_structured_grid(double x0, double y0, double z0, int nx, int ny, int nz, double dx);
    //int visualization_save(std::string& filename);

    int simplify(Polyhedron& boundary);
    int run(CGALSettings& settings);

    struct activity_sphere {
        K::Point_3* centre;
        float r;
        int i;
    };
}

#endif
