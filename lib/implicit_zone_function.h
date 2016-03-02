#ifndef IMPLICIT_ZONE_FUNCTION_H
#define IMPLICIT_ZONE_FUNCTION_H

#include "mesher_cgal.h"
#include "zone_priority_sorting.h"

#include <CGAL/Point_inside_polyhedron_3.h>
typedef CGAL::Point_inside_polyhedron_3<Polyhedron,K> Point_inside_polyhedron;

typedef std::map< int, Polyhedron* > region_ip_map;
typedef std::map< int, Point_inside_polyhedron*, ZonePrioritySorting > zone_pip_map;
typedef std::map< int, struct mesherCGAL::activity_sphere* > zone_activity_sphere_map;

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

template<typename K >
class Implicit_zone_function
{
public:
    typedef int                    return_type;
    typedef typename K::Point_3    Point_3;
    typedef typename K::FT         FT;
    typedef typename std::vector< std::pair<int, Point_inside_polyhedron*> > pip_vector;

    Implicit_zone_function(
            region_ip_map& region_ips,
            const Point_3* centre,
            float radius,
            int* organ_id,
            int* extent_id,
            std::vector<int>& vessels,
            std::vector<int>& needles,
            zone_pip_map& zone_pips,
            zone_activity_sphere_map& zone_activity_spheres,
            int default_zone=1
    ) : zone_pips_(zone_pips), default_zone_(default_zone), organ_id_(organ_id), extent_id_(extent_id),
        zone_activity_spheres_(zone_activity_spheres)
    {
      Polyhedron* organ = NULL;
      organ_pip_ = NULL;
      radius_ = radius;
      centre_ = centre;

      extent_tree_ = new Tree();

      if (organ_id != NULL) {
        organ = new Polyhedron(*region_ips[*organ_id]);
        extent_tree_->insert(organ->facets_begin(), organ->facets_end(), *organ);
        organ_pip_= new Point_inside_polyhedron(*organ);
      }

      extent_tree_->accelerate_distance_queries();

      BOOST_FOREACH(std::vector<int>::value_type& id, vessels) {
          Tree* tree = new Tree(region_ips[id]->facets_begin(), region_ips[id]->facets_end(), *region_ips[id]);
          tree->accelerate_distance_queries();
          trees_.push_back(tree);
          vessels_pip_.push_back(std::pair<int, Point_inside_polyhedron*>(id, new Point_inside_polyhedron(*region_ips[id])));
      }

      BOOST_FOREACH(std::vector<int>::value_type& id, needles) {
          Tree* tree = new Tree(region_ips[id]->facets_begin(), region_ips[id]->facets_end(), *region_ips[id]);
          tree->accelerate_distance_queries();
          trees_.push_back(tree);
          needles_pip_.push_back(std::pair<int, Point_inside_polyhedron*>(id, new Point_inside_polyhedron(*region_ips[id])));
      }
    }

    ~Implicit_zone_function()
    {
    }

    return_type signed_domain(const Point_3& p, const bool = true) const
    {
        int zone = default_zone_;

        if (CGAL::squared_distance(*centre_, p) > radius_ * radius_)
            return *extent_id_ == 0 ? 0 : -*extent_id_;

        if (organ_pip_ != NULL) {
            if (!(*organ_pip_)(p))
                return *organ_id_ == 0 ? -1000 : -*organ_id_; //FIXME: replace with something sensible (for when default_zone==0)
            else
                zone = *organ_id_;
        }

        BOOST_FOREACH(const pip_vector::value_type& pip, vessels_pip_) {
            if ((*pip.second)(p))
                return -pip.first;
        }
        BOOST_FOREACH(const pip_vector::value_type& pip, needles_pip_) {
            if ((*pip.second)(p))
                return -pip.first;
        }

        BOOST_FOREACH(const zone_pip_map::value_type& zone_pair, zone_pips_) {
            if ((*zone_pair.second)(p)) {
                zone = zone_pair.first;
                zone_activity_sphere_map::iterator activity_sphere = zone_activity_spheres_.find(zone_pair.first);
                if (activity_sphere != zone_activity_spheres_.end())
                {
                    Point_3* centre = activity_sphere->second->centre;
                    float radius = activity_sphere->second->r;
                    if (CGAL::squared_distance(*centre, p) > radius * radius)
                    {
                        return -activity_sphere->second->i;
                    }
                }
                break;
            }
        }
        if (default_zone_ == 0)
            zone += 1;

        return zone;
    }

    return_type operator()(const Point_3& p, const bool = true) const
    {
        return_type i = signed_domain(p);
        return i;
        //return i < 0 ? 0 : i;
    }

protected:
    zone_pip_map& zone_pips_;
    int default_zone_, *organ_id_, *extent_id_;

    const Point_3* centre_;
    float radius_;
    Tree *extent_tree_;
    std::vector<Tree*> trees_;
    Point_inside_polyhedron *extent_pip_, *organ_pip_;
    pip_vector vessels_pip_, needles_pip_;
    zone_activity_sphere_map& zone_activity_spheres_;
};

#endif
