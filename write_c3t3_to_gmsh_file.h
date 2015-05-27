/**
 * mesher_cgal
 *
 * Copyright (C) 2013-  NUMA Engineering Ltd. (see AUTHORs file)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Copyright (c) 2008-2009  GeometryFactory (France).
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 * You can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 * $URL: svn+ssh://sloriot@scm.gforge.inria.fr/svn/cgal/branches/next/Mesh_3/include/CGAL/IO/Complex_3_in_triangulation_3_to_vtk.h $
 * $Id: Complex_3_in_triangulation_3_to_vtk.h 67117 2012-01-13 18:14:48Z lrineau $
 *
 *
 * Author(s)     : Laurent Rineau
 */
// ACKNOWLEDGEMENT: This approach was inspired by David Bernstein's post to cgal-discuss

#ifndef Mesh_3_example_write_c3t3_to_gmsh_file_h
#define Mesh_3_example_write_c3t3_to_gmsh_file_h

namespace CGAL {
    template<class C3t3, class Polyhedron>
    bool write_c3t3_to_gmsh_file(const C3t3 &c3t3, std::map<typename C3t3::Triangulation::Facet,int>& boundary_indices, const std::string &file_name, bool decrement_zone=false, bool number_from_zero=false)
    {
        typedef typename C3t3::Triangulation Tr;
        typedef typename C3t3::Cells_in_complex_iterator Cell_iterator;
        typedef typename C3t3::Facets_in_complex_iterator Facet_iterator;
        typedef typename Tr::Finite_vertices_iterator Vertex_iterator;
        typedef typename Polyhedron::Halfedge_around_facet_circulator Halfedge_circulator;

        // Domain
        typedef Exact_predicates_inexact_constructions_kernel K;
        typedef K::FT FT;
        typedef K::Point_3 Point;

        // check that file extension is "msh"
        CGAL_assertion(file_name.substr(file_name.length()-4,4) == ".msh");

        // open file
        std::ofstream gmsh_file(file_name.c_str());

        // header
        gmsh_file << "$MeshFormat" << std::endl << "2.0 0 8" << std::endl << "$EndMeshFormat" << std::endl;

        // write mesh
        Tr t = c3t3.triangulation();
        int num_vertices = t.number_of_vertices();
        int num_facets = c3t3.number_of_facets_in_complex();
        int num_cells = c3t3.number_of_cells_in_complex();
        std::cout << num_facets << " - " << num_cells << std::endl;

        Cell_iterator it;

        // Write vertices
        std::map<Point, int> V;
        int i=0, j=0;
        int i_delta = number_from_zero ? 0 : 1;
        for (it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it)
            for (j = 0 ; j < 4 ; j++)
                    if (V.count(it->vertex(j)->point()) == 0) {
                            i++;
                            V[it->vertex(j)->point()] = -1;
                    }

        num_vertices = i;
        gmsh_file << "$Nodes" << std::endl;
        gmsh_file << num_vertices << std::endl;

	i = 0;
        for (Vertex_iterator it=t.finite_vertices_begin(); it != t.finite_vertices_end(); ++it)
        {
	    if (V.count(it->point()) == 0) {
		    std::cout << " (removed redundant vertex: " << i << ")" << std::endl;
		    continue;
	    }
	    gmsh_file << i + i_delta << " " << it->point().x() << " " << it->point().y() << " " << it->point().z() << std::endl;
            V[it->point()] = i;
            ++i;
        }

        gmsh_file << "$EndNodes" << std::endl;

        // Write tetrahedra
        gmsh_file << "$Elements" << std::endl;

        gmsh_file << num_cells + num_facets;
        gmsh_file << std::endl;

	std::vector<int> indices(3, 0);

        int elt = number_from_zero ? -1 : 0;
        Facet_iterator fit;
        for (fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit)
        {
          elt++;
          int j=-1;

          for (int i = 0; i < 4; ++i)
            if (i != fit->second)
                indices[++j] = V[(*fit).first->vertex(i)->point()];
          //if ( fit->second%2 == 1 ) std::swap(indices[0],indices[1]); RMV
          int subId = boundary_indices[(*fit)];
	  std::sort(indices.begin(), indices.end());
          gmsh_file << elt << " 2 2 " << subId << " " << subId << " " << indices[0] + i_delta <<" " << indices[1] + i_delta <<" " << indices[2] + i_delta << std::endl;
        }

        for (it = c3t3.cells_in_complex_begin(); it != c3t3.cells_in_complex_end(); ++it)
        {
            int subId = c3t3.subdomain_index(it);
            if (decrement_zone)
                 subId--;
            gmsh_file << ++elt << " 4 2 " << subId << " " << subId << " ";
            gmsh_file << V[it->vertex(0)->point()] + i_delta << " ";
            gmsh_file << V[it->vertex(1)->point()] + i_delta << " ";
            gmsh_file << V[it->vertex(2)->point()] + i_delta << " ";
            gmsh_file << V[it->vertex(3)->point()] + i_delta << std::endl;
        }

        gmsh_file << "$EndElements" << std::endl;

        return true;
    }
} // end namespace CGAL

#endif
