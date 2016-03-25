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

#include "Signed_mesh_domain_3.h"
#include "implicit_zone_function.h"
#include "proximity_domain_3.h"

namespace mesherCGAL {
    // MesherCGAL members
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
        std::vector< int > polyhedral_zones;
        std::vector< std::unique_ptr<Zone> > _zones;

        region_ep_map _region_eps;
        FT _bbox_radius;
        Tree *_boundary_tree;
        std::unique_ptr<Iso_cuboid> _bbox_p;
        Proximity_domain_field_3<Mesh_domain::R,Mesh_domain::Index>* _pdf;
        Mesh_domain_implicit* _domain;

        public:
            MesherCGAL(const CGALSettings& settings) :
                _settings(settings), _centre(NULL),
                _boundary_tree(NULL), _pdf(NULL),
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
