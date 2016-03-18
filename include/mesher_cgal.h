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

// tree
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Side_of_triangle_mesh.h>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

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

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef CGAL::Side_of_triangle_mesh<Polyhedron,K> Side_of_triangle_mesh;

namespace mesherCGAL {

    //int visualization_set_allocation_order(int ix, int order, double cl, double needle_dist, double x, double y, double z);
    //int visualization_create_structured_grid(double x0, double y0, double z0, int nx, int ny, int nz, double dx);
    //int visualization_save(std::string& filename);
    typedef std::map< int, std::shared_ptr<Polyhedron> > region_ip_map;
    typedef std::map< int, std::shared_ptr<Exact_polyhedron> > region_ep_map;

    int simplify(Polyhedron& boundary);

    struct activity_sphere {
        K::Point_3* centre;
        float r;
        int i;
    };

    class Zone {
        private:
            /* If this ever ends up in a scenario where zone is destroyed, switch to unique_ptr and
             * check performance impact */
            struct activity_sphere* activity_sphere;

        protected:
            int _id;
            float _cl;
            float _priority;

        public:
            Zone(int id, float cl, float priority) :
                activity_sphere(NULL), _id(id), _cl(cl), _priority(priority) {}
            virtual ~Zone() {}
            bool set_activity_sphere(float x, float y, float z, float r, int i);
            void print_activity_sphere();
            int get_id() const { return _id; }
            float get_priority() const { return _priority; }
            float get_cl() const { return _cl; }
            virtual bool is_container() const = 0;
            virtual std::shared_ptr<Exact_polyhedron> exact_polyhedron() = 0;

            int contains(const K::Point_3& p, bool check_activity_sphere=true) const {
                if (!is_container())
                    return 0;

                int id = contains_all(p) ? 0 : _id;

                if (id && check_activity_sphere) {
                    int inactivity_index;
                    if ((inactivity_index = outside_activity_sphere(p)) != 0)
                        return inactivity_index;
                }

                return id;
            }
            virtual bool has_tree() const = 0;
            virtual float squared_distance(const K::Point_3& p) const = 0;
            virtual bool add_tree() = 0;
            virtual bool load(std::string filename) = 0;

        protected:
            /* Ignores inactivity */
            virtual bool contains_all(const K::Point_3& p) const = 0;
            int outside_activity_sphere(const K::Point_3& p) const {
                if (activity_sphere) {
                    Point* centre = activity_sphere->centre;
                    float radius = activity_sphere->r;
                    if (CGAL::squared_distance(*centre, p) > radius * radius)
                    {
                        return -activity_sphere->i;
                    }
                }

                return 0;
            }
    };

    class PolyhedralZone : public Zone {
        Tree* _tree;
        std::shared_ptr<Side_of_triangle_mesh> _pip;
        std::shared_ptr<Polyhedron> _ip;
        std::shared_ptr<Exact_polyhedron> _ep;

        public:
            PolyhedralZone(int id, float cl, float priority) : Zone(id, cl, priority),
                _tree(NULL), _ip(NULL), _ep(NULL) {}
            bool is_container() const { return (bool)_pip; }
            std::shared_ptr<Exact_polyhedron> exact_polyhedron() { return _ep; }

            bool has_tree() const { return _tree; }
            float squared_distance(const K::Point_3& p) const {
                if (!has_tree())
                    return -1.0;
                return _tree->squared_distance(p);
            }
            bool add_tree() {
                if (!_ip || has_tree())
                    return false;

                _tree = new Tree();
                _tree->insert(_ip->facets_begin(), _ip->facets_end(), *_ip);

                return true;
            }
            bool load(std::string filename);

        protected:
            bool contains_all(const K::Point_3& p) const {
                return (*_pip)(p) == CGAL::ON_UNBOUNDED_SIDE;
            }

        private:
            bool create_polyhedra(std::shared_ptr<Exact_polyhedron> ep);
    };

    class UnstructuredGridZone : public Zone {
        vtkSmartPointer<vtkUnstructuredGrid> _grid;

        public:
            UnstructuredGridZone(int id, float cl, float priority) : Zone(id, cl, priority) {}
            bool is_container() const { return true; };
            std::shared_ptr<Exact_polyhedron> exact_polyhedron() { return NULL; };

            bool has_tree() const { return false; };
            float squared_distance(const K::Point_3&) const { return 0; };
            bool add_tree() { return false; };
            bool load(std::string filename);

        protected:
            bool contains_all(const K::Point_3& p) const;
    };
}

#endif
