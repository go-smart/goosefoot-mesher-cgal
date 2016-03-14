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

#ifndef MESHER_CGAL_APP_H
#define MESHER_CGAL_APP_H

#include "mesher_cgal.h"

#include "zone_priority_sorting.h"
#include "Signed_mesh_domain_3.h"
#include "implicit_zone_function.h"
#include "proximity_domain_3.h"

namespace mesherCGAL {
    // MesherCGAL members
    typedef std::map< int, Polyhedron* > region_ip_map;
    typedef std::map< int, Exact_polyhedron* > region_ep_map;
    typedef std::map< int, Polyhedron* > zone_ip_map;
    typedef std::map< int, struct activity_sphere* > zone_activity_sphere_map;
    typedef CGAL::Signed_mesh_domain_3< Implicit_zone_function<K>, K > Mesh_domain_implicit;
    typedef K::Iso_cuboid_3 Iso_cuboid;

    class MesherCGAL {
        CGALSettings _settings;
        Point* _centre;
        std::map<typename C3t3::Triangulation::Facet,int> _boundary_indices;
        region_ip_map _region_ips;
        C3t3 _c3t3;
        std::vector< int > vessels;
        std::vector< int > needles;
        std::vector< int > zones, polyhedral_zones;
        zone_cls_map _zone_priorities;
        zone_cls_map _zone_cls;
        zone_pip_map _zone_pips;
        region_ep_map _region_eps;
        zone_ip_map _zone_ips;
        zone_activity_sphere_map _zone_activity_spheres;
        FT _bbox_radius;
        Tree *_boundary_tree;
        Iso_cuboid* _bbox_p;
        Proximity_domain_field_3<Mesh_domain::R,Mesh_domain::Index>* _pdf;
        Mesh_domain_implicit* _domain;

        public:
            MesherCGAL(const CGALSettings& settings) :
                _settings(settings), _centre(NULL),
                _zone_pips(ZonePrioritySorting(_zone_priorities)),
                _boundary_tree(NULL), _bbox_p(NULL), _pdf(NULL),
                _domain(NULL)
            {}

            int init();
            int setup_regions();
            int calculate_bbox();
            int setup_domain_field();
            int mesh();
            int label_boundaries();
            int output();

            int run();
    };
}

#endif
