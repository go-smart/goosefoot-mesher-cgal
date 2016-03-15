#ifndef IMPLICIT_ZONE_FUNCTION_H
#define IMPLICIT_ZONE_FUNCTION_H

#include "mesher_cgal.h"

#include <CGAL/Point_inside_polyhedron_3.h>
typedef CGAL::Point_inside_polyhedron_3<Polyhedron,K> Point_inside_polyhedron;

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

namespace mesherCGAL {
    typedef std::vector< Zone > zones_vec;

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
                bool use_organ,
                int organ_id,
                bool use_extent,
                int extent_id,
                std::vector<int>& vessels,
                std::vector<int>& needles,
                zones_vec& zones,
                int default_zone=1
        ) : zones_(zones), default_zone_(default_zone), use_organ_(use_organ), organ_id_(organ_id), use_extent_(use_extent), extent_id_(extent_id)
        {
          Polyhedron* organ = NULL;
          organ_pip_ = NULL;
          radius_ = radius;
          centre_ = centre;

          extent_tree_ = new Tree();

          if (use_organ) {
            organ = new Polyhedron(*region_ips[organ_id]);
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
                return extent_id_ == 0 ? 0 : -extent_id_;

            if (organ_pip_ != NULL) {
                if (!(*organ_pip_)(p))
                    return organ_id_ == 0 ? -1000 : -organ_id_; //FIXME: replace with something sensible (for when default_zone==0)
                else
                    zone = organ_id_;
            }

            BOOST_FOREACH(const pip_vector::value_type& pip, vessels_pip_) {
                if ((*pip.second)(p))
                    return -pip.first;
            }
            BOOST_FOREACH(const pip_vector::value_type& pip, needles_pip_) {
                if ((*pip.second)(p))
                    return -pip.first;
            }

            for (auto&& z : zones_) {
                if (z.contains(p) != 0) {
                    zone = z.get_id();
                    if (zone > 0 && default_zone_ == 0)
                        zone += 1;
                    return zone;
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
        zones_vec& zones_;
        int default_zone_;
        bool use_organ_;
        int organ_id_;
        bool use_extent_;
        int extent_id_;

        const Point_3* centre_;
        float radius_;
        Tree *extent_tree_;
        std::vector<Tree*> trees_;
        Point_inside_polyhedron *extent_pip_, *organ_pip_;
        pip_vector vessels_pip_, needles_pip_;
    };
}

#endif
