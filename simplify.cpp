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

/*
 * This file is adapted under the very permissive CGAL-4.3 examples licence
 * Note that there were therefore unnamed upstream authors. The license
 * is reproduced as follows:
 *     Copyright (c) 1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007
 *     Utrecht University (The Netherlands),
 *     ETH Zurich (Switzerland),
 *     INRIA Sophia-Antipolis (France),
 *     Max-Planck-Institute Saarbruecken (Germany),
 *     and Tel-Aviv University (Israel).  All rights reserved.
 *
 *     Permission is hereby granted, free of charge, to any person obtaining a
 *     copy of this software and associated documentation files (the
 *     "Software"), to deal in the Software without restriction, including
 *     without limitation the rights to use, copy, modify, merge, publish,
 *     distribute, sublicense, and/or sell copies of the Software, and to
 *     permit persons to whom the Software is furnished to do so, subject to
 *     the following conditions:
 *
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 *     OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 *     CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 *     TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * */

#include <boost/variant.hpp>
#include <boost/tuple/tuple.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
namespace SMS = CGAL::Surface_mesh_simplification;

int simplify(Polyhedron& boundary) {
	SMS::Count_stop_predicate<Polyhedron> stop(1000);
#if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(4,5,0)
	int r = SMS::edge_collapse
		(boundary
		, stop
		, CGAL::vertex_index_map(get(CGAL::vertex_external_index, boundary))
		.halfedge_index_map(get(CGAL::halfedge_external_index, boundary))
		.get_cost(SMS::Edge_length_cost <Polyhedron>())
		.get_placement(SMS::Midpoint_placement<Polyhedron>())
		);
#else
	int r = SMS::edge_collapse(boundary, stop,
		CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index, boundary))
		.edge_index_map(boost::get(CGAL::edge_external_index, boundary)));
#endif

	return r;
}

