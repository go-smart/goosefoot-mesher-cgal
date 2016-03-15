#ifndef ZONE_PRIORITY_SORTING_H
#define ZONE_PRIORITY_SORTING_H

#include "mesher_cgal.h"

#include <CGAL/Point_inside_polyhedron_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef std::map< int, float > zone_cls_map, zone_priorities_map;
typedef std::map< int, Tree* > zone_tree_map;

typedef CGAL::Point_inside_polyhedron_3<Polyhedron,K> Point_inside_polyhedron;

namespace mesherCGAL {
    class ZonePrioritySorting
    {
        public:
            ZonePrioritySorting() {}
            bool operator() (const Zone& z, const Zone& w) const
            { return fabs(z.get_priority() - w.get_priority()) < 1e-9 ? (z.get_id() < w.get_id()) : z.get_priority() < w.get_priority(); }
    };

    typedef std::map< int, Point_inside_polyhedron*, ZonePrioritySorting > zone_pip_map;
}

#endif
