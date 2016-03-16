#ifndef ZONE_PRIORITY_SORTING_H
#define ZONE_PRIORITY_SORTING_H

#include "mesher_cgal.h"

#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef CGAL::Side_of_triangle_mesh<Polyhedron,K> Side_of_triangle_mesh;

namespace mesherCGAL {
    class ZonePrioritySorting
    {
        public:
            ZonePrioritySorting() {}
            bool operator() (const std::unique_ptr<Zone>& z, const std::unique_ptr<Zone>& w) const
            { return fabs(z->get_priority() - w->get_priority()) < 1e-9 ? (z->get_id() < w->get_id()) : z->get_priority() < w->get_priority(); }
    };
}

#endif
