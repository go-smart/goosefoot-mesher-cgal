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

class ZonePrioritySorting
{
    private:
        zone_cls_map& _zone_priorities;
    public:
        ZonePrioritySorting(zone_cls_map& zone_priorities) : _zone_priorities(zone_priorities) {}
        bool operator() (const int& i, const int& j) const
        { return fabs(_zone_priorities[i] -_zone_priorities[j]) < 1e-9 ? (i < j) : _zone_priorities[i] < _zone_priorities[j]; }
};

typedef std::map< int, Point_inside_polyhedron*, ZonePrioritySorting > zone_pip_map;

#endif
