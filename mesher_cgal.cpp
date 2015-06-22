/**
 * mesher_cgal
 *
 * Copyright (C) 2013-  NUMA Engineering Ltd. (see AUTHORS file)
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

#include <sstream>
#include <functional>
#include <sys/stat.h>

#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Point_inside_polyhedron_3.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include "Signed_mesh_domain_3.h"
// IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/utils_classes.h>

// surface mesh
#include <CGAL/Polyhedron_3.h>

// nef
#include <CGAL/Nef_polyhedron_3.h> 
#include <CGAL/Nef_3/SNC_indexed_items.h>

#include <CGAL/bounding_box.h>

#include "mesher_cgal.h"
#include "mesher_cgal.int.h"
#include "write_c3t3_to_gmsh_file.h"

#include <google/protobuf/text_format.h>

double EPS = 0.00001;

using namespace mesherCGAL;

template <typename Refs>
struct MyFace : public CGAL::HalfedgeDS_face_base<Refs> {
    int boundary_index;
};

struct MyItems : public CGAL::Polyhedron_items_3 {
    template <typename Refs, typename Traits>
    struct Face_wrapper {
        typedef MyFace<Refs> Face;
    };
};

// Tree
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
typedef K::Triangle_3 Triangle;
typedef K::Plane_3 Plane;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef K::Iso_cuboid_3 Iso_cuboid;
typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef CGAL::AABB_polyhedron_triangle_primitive<Exact_Kernel,Exact_polyhedron> Exact_Primitive;
typedef CGAL::AABB_traits<Exact_Kernel, Exact_Primitive> Exact_Traits;
typedef CGAL::AABB_tree<Exact_Traits> Exact_Tree;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

//#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
//#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

//namespace SMS = CGAL::Surface_mesh_simplification;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

struct activity_sphere {
    K::Point_3* centre;
    float r;
    int i;
};

typedef std::map< int, struct activity_sphere* > zone_activity_sphere_map;
typedef std::map< int, int > region_boundary_count_map;
typedef std::map< int, std::string > region_string_map;
typedef std::map< int, Exact_polyhedron* > region_ep_map;
typedef std::map< int, Polyhedron* > region_ip_map;
typedef std::map< int, Exact_Tree* > region_tree_map;

typedef std::map< int, std::string > zone_string_map;
typedef std::map< int, Polyhedron* > zone_ip_map;
typedef std::map< int, Exact_polyhedron* > zone_ep_map;
typedef std::map< int, float > zone_cls_map, zone_priorities_map;
typedef std::map< int, Tree* > zone_tree_map;

#include <map>
#include <vector>
#include <utility>

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

template <typename Gt, typename Index_>
class Proximity_domain_field_3
{
public:
  typedef typename Gt::FT         FT;
  typedef typename Gt::Point_3    Point_3;
  typedef Index_                  Index;
  
private:
  // Map to store field values
  int minx, miny, minz, nx, ny, nz;
  FT length_scale;
  CGAL::Max<FT> max_;
  CGAL::Min<FT> min_;
  
public:
  /// Constructor
  Proximity_domain_field_3(Tree* tree, double near_cl, double far_cl, zone_cls_map& zone_cls, float zone_radius, float centre_radius,
                const Point* centre, Tree* needle_tree, zone_tree_map& zone_trees, float granularity, const Iso_cuboid& bbox,
                zone_pip_map* zone_pips) :
                        tree_(tree), near_cl_(near_cl), far_cl_(far_cl), zone_cls_(zone_cls), zone_radius_(zone_radius),
			centre_radius_(centre_radius),
                        centre_(centre), needle_tree_(needle_tree), zone_trees_(zone_trees), granularity_(granularity),
                        bbox_(bbox), checks_(0), total_checks_(0), values_size_(0), zone_pips_(zone_pips),
                        max_(CGAL::Max<FT>()), min_(CGAL::Min<FT>())
    {
        max_dist[0] = 0.; max_dist[1] = 0.; max_dist[2] = 0.;

        minx = (int)(bbox.xmin()/granularity_); nx = 1 + (int)(bbox.xmax()/granularity_) - minx;
        miny = (int)(bbox.ymin()/granularity_); ny = 1 + (int)(bbox.ymax()/granularity_) - miny;
        minz = (int)(bbox.zmin()/granularity_); nz = 1 + (int)(bbox.zmax()/granularity_) - minz;
        length_scale = std::pow(bbox.volume(), 1./3);
        
        //visualization_create_structured_grid(minx, miny, minz, nx, ny, nz, granularity_);

        values_ = (FT*)calloc(sizeof(FT), nx*ny*nz);
	max_cl_ = -1.;
    }
  //RMV : insert destructor to clear values_ allocations

  FT to_cl(FT closest_cl, FT dist) const
  {
      FT scaling = 1;

      if (dist <= 4 * closest_cl * closest_cl + 1e-10)
          scaling = 0;
      else
          scaling *= pow((dist - closest_cl * closest_cl * 4)/(length_scale * length_scale * 0.25), 0.5);

      return (far_cl_ - (far_cl_-closest_cl)*exp(-scaling));;
  }

  /// Returns size
  FT operator()(const Point_3& p, const int dim, const Index& index) const
  {
      std::vector<int> *tup = new std::vector<int>;
      int x = (int)(p.x()/granularity_);
      int y = (int)(p.y()/granularity_);
      int z = (int)(p.z()/granularity_);
      int ix = nx*(ny*(min_(max_(z - minz, 0), nz)) + min_(max_(y - miny, 0), ny)) + min_(max_(x - minx, 0), nx);

      FT cl = values_[ix];

      checks_++;

      if ( checks_ % 10000 == 0 ) {
          total_checks_ += checks_;
          std::cout << "Cached " << values_size_ << " blocks (miss rate " << values_size_*100./total_checks_ << "%) - max. cl: " << max_cl_ << " - max. distances from key points/surfaces: " << max_dist[0] << " " << max_dist[1] << " " << max_dist[2] << std::endl;
          checks_ = 0;

          //if (values_.size() > 0) {
          //    typename Values::const_iterator i = values_.begin();
          //    char filename[200];
          //    sprintf(filename, "/tmp/points.%06d.txt", total_checks_);
          //    FILE* f = fopen(filename, "w");
          //    fprintf(f, "x,y,z\n");
          //    for ( ; i != values_.end() ; ++i ) {
          //        const std::vector<int> &curtup = (i->first);
          //        fprintf(f, "%d, %d, %d\n", curtup[0], curtup[1], curtup[2]);
          //    }
          //    fclose(f);
          //}
      }

      if ( cl > 1e-20 ) {
          return cl;
      }

      FT scaling = 1., dist = -1.;
      FT closest_cl = near_cl_;

      if (tree_ != NULL) {
          FT tree_dist = tree_->squared_distance(p);
          if (dist < 0 || tree_dist < dist) {
              dist = tree_dist;
              closest_cl = near_cl_;
          }
          //max_dist[0] = CGAL::max(max_dist[0], tree_dist);
      }

      FT needle_dist = 0.;
      if (needle_tree_ != NULL) {
          needle_dist = needle_tree_->squared_distance(p);
          //needle_dist = CGAL::min(needle_dist - near_cl_ * 2, 0.0); /* Introduce near_cl region around needle */
          if (dist < 0 || needle_dist < dist) {
              dist = needle_dist;
              closest_cl = near_cl_;
          }

      }

      BOOST_FOREACH(const zone_tree_map::value_type& zone_pair, zone_trees_) {
          FT zone_dist = -1;
          Tree* zone_tree = zone_pair.second;
          if (zone_pips_ != NULL && (*(*zone_pips_)[zone_pair.first])(p)) {
              zone_dist = 0;
          }
          else {
              zone_dist = zone_tree->squared_distance(p);
              if (zone_dist > zone_radius_ * zone_radius_) 
                  zone_dist = zone_dist + zone_radius_ * zone_radius_ - 2. * sqrt(zone_dist) * zone_radius_;
              else zone_dist = 0.;
              if (zone_dist < 0) zone_dist = 0.;
          }
          //zone_dist = CGAL::min(zone_dist - near_cl_ * 2, 0.0); /* Introduce near_cl region around zone */
          if (dist < 0 || (to_cl(closest_cl, dist) > to_cl(zone_cls_[zone_pair.first], zone_dist))) {
              dist = zone_dist;
              closest_cl = zone_cls_[zone_pair.first];
          }
      }

      if (centre_ != NULL) {
          FT centre_dist = squared_distance(*centre_, p);
          if (centre_dist > centre_radius_ * centre_radius_)
              centre_dist = centre_dist + centre_radius_ * centre_radius_ - 2. * sqrt(centre_dist) * centre_radius_;
          else centre_dist = 0.;
          if (centre_dist < 0) centre_dist = 0.;
          if (dist < 0 || centre_dist < dist) {
              dist = centre_dist;
              closest_cl = near_cl_;
          }
      }

      if (dist < 0) {
          cl = far_cl_;
      } else {
          //if (dist <= 4 * closest_cl * closest_cl + 1e-10)
          //    scaling = 0;
          //else
          //    scaling *= pow((dist - closest_cl * closest_cl * 4)/(length_scale * length_scale * 0.25), 0.5);
          cl = CGAL::min(far_cl_, to_cl(closest_cl, dist));
	  if (cl > max_cl_) max_cl_ = cl;
      }

      //if (zone_pips_ != NULL) {
      //    int zone = -1;
      //    BOOST_FOREACH(const zone_pip_map::value_type& zone_pair, *zone_pips_) {
      //        if ((*zone_pair.second)(p)) {
      //            zone = zone_pair.first;
      //            break;
      //        }
      //    }
      //    if (zone > -1 && zone_cls_[zone] < cl) {
      //      cl = zone_cls_[zone];
      //      if (zone_cls_[zone] < 0.01)
      //          std::cout << "*" << std::endl;
      //    }
      //}

      //FT scaling = CGAL::min(10*needle_dist/length_scale, 4.) * (1 / (.1 + dist/length_scale));
      //printf("[%lf %lf %lf]\n", needle_dist, dist, length_scale);

      //cl = CGAL::min(far_cl_, (far_cl_ - (far_cl_-near_cl_)*exp(-scaling)));
      values_size_++;

      values_[ix] = cl;
      //visualization_set_allocation_order(ix, values_size_, cl, needle_dist, p.x(), p.y(), p.z());

      return cl;
  }

  mutable float max_cl_;//RMV
private:
  Tree* tree_;
  mutable FT* values_;
  FT near_cl_, far_cl_;
  zone_cls_map& zone_cls_;
  FT zone_radius_, centre_radius_;
  const Point* centre_;
  Tree* needle_tree_;
  zone_tree_map zone_trees_;
  FT granularity_;
  const Iso_cuboid& bbox_;
  mutable int checks_, total_checks_, values_size_;
  mutable FT max_dist[3];
  zone_pip_map* zone_pips_;
};

// See : http://cgal-discuss.949826.n4.nabble.com/Example-Convert-polyhedron-from-one-kernel-to-another-td4514497.html
// Can be used to convert polyhedron from exact to inexact and vice-versa
template <class Polyhedron_input,
class Polyhedron_output>
struct Copy_polyhedron_to
        : public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
        Copy_polyhedron_to(const Polyhedron_input& in_poly)
                : in_poly(in_poly) {}

        void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
        {
                typedef typename Polyhedron_output::HalfedgeDS Output_HDS;
                typedef typename Polyhedron_input::HalfedgeDS Input_HDS;

                CGAL::Polyhedron_incremental_builder_3<Output_HDS> builder(out_hds);

                typedef typename Polyhedron_input::Vertex_const_iterator Vertex_const_iterator;
                typedef typename Polyhedron_input::Facet_const_iterator  Facet_const_iterator;
                typedef typename Polyhedron_input::Halfedge_around_facet_const_circulator HFCC;

                builder.begin_surface(in_poly.size_of_vertices(),
                        in_poly.size_of_facets(),
                        in_poly.size_of_halfedges());

                for(Vertex_const_iterator
                        vi = in_poly.vertices_begin(), end = in_poly.vertices_end();
                        vi != end ; ++vi)
                {
                        typename Polyhedron_output::Point_3 p(::CGAL::to_double( vi->point().x()),
                                ::CGAL::to_double( vi->point().y()),
                                ::CGAL::to_double( vi->point().z()));
                        builder.add_vertex(p);
                }

                typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
                Index index( in_poly.vertices_begin(), in_poly.vertices_end());

                for(Facet_const_iterator
                        fi = in_poly.facets_begin(), end = in_poly.facets_end();
                        fi != end; ++fi)
                {
                        HFCC hc = fi->facet_begin();
                        HFCC hc_end = hc;
                        builder.begin_facet ();
                        do {
                                builder.add_vertex_to_facet(index[hc->vertex()]);
                                ++hc;
                        } while( hc != hc_end);
                        builder.end_facet();
                }
                builder.end_surface();
        } // end operator()(..)
private:
        const Polyhedron_input& in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_B, class Poly_A>
void poly_copy(Poly_B& poly_b, const Poly_A& poly_a)
{
        poly_b.clear();
        Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
        poly_b.delegate(modifier);
} 

int parse_region(region_string_map& region_files, const std::string& input_string) {
  std::vector<std::string> pair;
  boost::split(pair, input_string, boost::is_any_of(":"));

  int id = boost::lexical_cast<int>(pair[1]);
  region_files[id] = pair[0];

  return id;
}

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
        extent_tree_->insert(organ->facets_begin(), organ->facets_end());
        organ_pip_= new Point_inside_polyhedron(*organ);
      }

      extent_tree_->accelerate_distance_queries();

      BOOST_FOREACH(std::vector<int>::value_type& id, vessels) {
          Tree* tree = new Tree(region_ips[id]->facets_begin(), region_ips[id]->facets_end());
          tree->accelerate_distance_queries();
          trees_.push_back(tree);
          vessels_pip_.push_back(std::pair<int, Point_inside_polyhedron*>(id, new Point_inside_polyhedron(*region_ips[id])));
      }

      BOOST_FOREACH(std::vector<int>::value_type& id, needles) {
          Tree* tree = new Tree(region_ips[id]->facets_begin(), region_ips[id]->facets_end());
          tree->accelerate_distance_queries();
          trees_.push_back(tree);
          needles_pip_.push_back(std::pair<int, Point_inside_polyhedron*>(id, new Point_inside_polyhedron(*region_ips[id])));
      }
    }

    ~Implicit_zone_function()
    {
    }

    return_type operator()(const Point_3& p, const bool = true) const
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


typedef CGAL::Signed_mesh_domain_3< Implicit_zone_function<K>, K > Mesh_domain_implicit;
typedef CGAL::Mesh_triangulation_3<Mesh_domain_implicit>::type Exact_Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Exact_Tr> Exact_C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Exact_Mesh_criteria;

int mesherCGAL::run(CGALSettings& settings) {

  std::string output_prefix("gsmesh.out");
  bool include_boundary_tree = false, omit_needle_tree = false, omit_zone_tree = false, include_needle = false,
       include_centre = false, solid_zone = false, mark_zone_boundaries = false;
  bool number_from_zero = false, tetrahedralize_only = false;
  float nearfield, farfield, zonefield;
  float granularity = 1., centre_radius = 1., bounding_radius = 50., zone_radius = 0.;
  region_string_map region_files;
  zone_string_map zone_files;
  region_ip_map region_ips;
  region_ep_map region_eps;
  zone_ip_map zone_ips;
  zone_ep_map zone_eps;
  zone_cls_map zone_cls;
  zone_cls_map zone_priorities;
  zone_activity_sphere_map zone_activity_spheres;
  std::vector< int > vessels;
  std::vector< int > needles;
  std::vector< int > zones, polyhedral_zones;
  int organ, extent;
  int default_zone_id;

  tetrahedralize_only = settings.tetrahedralize_only();
  include_needle = settings.needles_size() > 0;

  std::string settings_string;
  google::protobuf::TextFormat::PrintToString(settings, &settings_string);

  if (settings.dump_settings()) {
      std::cout << "LOADED SETTINGS FROM PROTOBUF" << std::endl
          << "+++++++++++++++++++++++++++++++" << std::endl
          << settings_string
          << "+++++++++++++++++++++++++++++++" << std::endl << std::endl;
  }

  std::vector<float> coords;
  Point* centre = NULL;

  if (settings.centre_size() > 0) {
      if (settings.centre_size() != 3) {
          std::cerr << "Need an \"x,y,z\" argument to define centre point - need 3 coords not " << coords.size() << std::endl;
          abort();
      }

      centre = new Point(settings.centre(0), settings.centre(1), settings.centre(2));

      std::cout << "Centre is set at (" << centre->x() << ", " << centre->y() << ", " << centre->z() << ")" << std::endl;
  }

  bool output_medit, output_gmsh, output_vtk, verbose;

  output_medit = settings.output_medit();
  number_from_zero = settings.number_from_zero();
  verbose = settings.verbose();
  output_gmsh = settings.output_gmsh();
  output_vtk = settings.output_vtk();
  include_centre = settings.dense_centre();
  omit_needle_tree = settings.omit_needle_tree();
  include_boundary_tree = settings.boundary_tree();
  omit_zone_tree = settings.omit_zone_tree();

  nearfield = settings.near_field();
  farfield = settings.far_field();
  zonefield = settings.zone_field();
  zone_radius = settings.zone_radius();

  if (settings.has_granularity())
      granularity = settings.granularity();
  else
      granularity = 1.;

  output_prefix = settings.output_prefix();

  centre_radius = settings.centre_radius();
  default_zone_id = settings.tissue_id();
  bounding_radius = settings.bounding_radius();
  tetrahedralize_only = settings.tetrahedralize_only();
  if (settings.has_mark_zone_boundaries())
      mark_zone_boundaries = settings.mark_zone_boundaries();
  if (settings.has_solid_zone())
      solid_zone = settings.solid_zone();

  std::cout << "Finding regions : " << std::endl;

  bool include_organ = settings.has_organ_file();

  if (include_organ) {
      organ = settings.organ_index();
      region_files[organ] = settings.organ_file();
      std::cout << " - Organ : " << organ << std::endl;
  }

  extent = settings.extent_index();
  //region_files[extent] = settings.extent_file();
  std::cout << " - Extent : " << extent << std::endl;

  include_needle = settings.needles_size() > 0;
  if (include_needle) {
      std::cout << " - Needles : " << std::endl << "   \\ ";
      for (int i = 0 ; i < settings.needles_size() ; i++) {
          int id = settings.needles_index(i);
          region_files[id] = settings.needles(i);
          needles.push_back(id);
          std::cout << id << " ";
      }
      std::cout << std::endl;
  }

  if (settings.vessels_size()) {
      std::cout << " - Vessels : " << std::endl << "   \\ ";
      for (int i = 0 ; i < settings.vessels_size() ; i++) {
          int id = settings.vessels_index(i);
          region_files[id] = settings.vessels(i);
          vessels.push_back(id);
          std::cout << id << " ";
      }
      std::cout << std::endl;
  }

  if (settings.zones_size() || settings.zone_size()) {
      std::cout << " - Zones : " << std::endl << "   \\ ";
      for (int i = 0 ; i < settings.zones_size() ; i++) {
          int id = settings.zones_index(i);
          zone_files[id] = settings.zones(i);
          zones.push_back(id);
          polyhedral_zones.push_back(id);
          zone_cls[id] = zonefield;
          zone_priorities[id] = 0.0;
          std::cout << id << " ";
      }
      for (int i = 0 ; i < settings.zone_size() ; i++) {
          int id = settings.zone(i).index();
          zone_files[id] = settings.zone(i).file();
          zones.push_back(id);
          polyhedral_zones.push_back(id);
          std::cout << id;

          if (settings.zone(i).has_characteristic_length()) {
              zone_cls[id] = settings.zone(i).characteristic_length();
              std::cout << "(" << zone_cls[id] << ")";
          }
          else {
              zone_cls[id] = zonefield;
          }

          if (settings.zone(i).has_priority()) {
              zone_priorities[id] = settings.zone(i).priority();
              if (fabs(zone_priorities[id]) > 1e-10)
                  std::cout << "<" << zone_priorities[id] << ">";
          }
          else {
              zone_priorities[id] = 0.0;
          }

          if (settings.zone(i).has_activity() && settings.zone(i).has_inactivity_index()) {
              zone_activity_spheres[id] = new struct activity_sphere;
              double
                  x = settings.zone(i).activity().x(),
                  y = settings.zone(i).activity().y(),
                  z = settings.zone(i).activity().z();
              if (centre != NULL)
              {
                  x += centre->x();
                  y += centre->y();
                  z += centre->z();
              }
              zone_activity_spheres[id]->centre = new Point(x, y, z);
              zone_activity_spheres[id]->r = settings.zone(i).activity().r();
              zone_activity_spheres[id]->i = settings.zone(i).inactivity_index();
              zones.push_back(zone_activity_spheres[id]->i);
              std::cout << "[" << zone_activity_spheres[id]->i << ":" << zone_activity_spheres[id]->r << "]";
          }
          std::cout << " ";
      }
      std::cout << std::endl;
  }

  ZonePrioritySorting zps(zone_priorities);
  zone_pip_map zone_pips(zps);

  if (settings.has_matching_tolerance()) {
      EPS = settings.matching_tolerance();
  }

  Polyhedron inexact_structures_polyhedron, inexact_combined_polyhedron, inexact_boundary_polyhedron;

  Exact_polyhedron exact_structures_polyhedron;
  Exact_polyhedron exact_boundary_polyhedron;
  Exact_polyhedron exact_combined_polyhedron;

  std::cout << "Loading (boundary) regions : " << std::endl << " - ";
  BOOST_FOREACH(const region_string_map::value_type& region_pair, region_files) {
      Exact_polyhedron* ep = new Exact_polyhedron();
      PolyhedronUtils::readSurfaceFile(region_pair.second, *ep);
      
      region_eps[region_pair.first] = ep;

      std::cout << region_pair.first << " (" << ep->size_of_facets() << " facets) " << std::flush;
  }
  std::cout << std::endl;

  std::cout << "Loading zones : " << std::endl << " - ";
  BOOST_FOREACH(const zone_string_map::value_type& zone_pair, zone_files) {
      Exact_polyhedron* ep = new Exact_polyhedron();
      PolyhedronUtils::readSurfaceFile(zone_pair.second, *ep);
      
      zone_eps[zone_pair.first] = ep;
      if (mark_zone_boundaries)
          region_eps[zone_pair.first] = ep;

      std::cout << zone_pair.first << " (" << ep->size_of_facets() << " facets) " << std::flush;
  }
  std::cout << std::endl;

  if (include_organ) {
      Polyhedron* inexact_organ_polyhedron = new Polyhedron();
      region_ips[organ] = inexact_organ_polyhedron;
      poly_copy<Polyhedron,Exact_polyhedron>(*inexact_organ_polyhedron, *region_eps[organ]);
  }

  Iso_cuboid* bbox_p = NULL;

  if (include_organ) {
      bbox_p = new Iso_cuboid(CGAL::bounding_box(region_ips[organ]->points_begin(), region_ips[organ]->points_end()));
      if (centre == NULL)
          centre = new Point(((*bbox_p)[0].x() + (*bbox_p)[7].x()) / 2,
                             ((*bbox_p)[0].y() + (*bbox_p)[7].y()) / 2,
                             ((*bbox_p)[0].z() + (*bbox_p)[7].z()) / 2);
  } else if (centre != NULL) {
      bbox_p = new Iso_cuboid(
          centre->x() - bounding_radius,
          centre->y() - bounding_radius,
          centre->z() - bounding_radius,
          centre->x() + bounding_radius,
          centre->y() + bounding_radius,
          centre->z() + bounding_radius
      );
  } else {
      std::cerr << "No centre defined and no organ to guess from" << std::endl;
      abort();
  }
  const Iso_cuboid& bbox(*bbox_p);

  BOOST_FOREACH(std::vector<int>::value_type& id, vessels) {
      Polyhedron* inexact_boundary_polyhedron = new Polyhedron();
      region_ips[id] = inexact_boundary_polyhedron;
      poly_copy<Polyhedron,Exact_polyhedron>(*region_ips[id], *region_eps[id]);
  }

  BOOST_FOREACH(std::vector<int>::value_type& id, needles) {
      Polyhedron* inexact_boundary_polyhedron = new Polyhedron();
      region_ips[id] = inexact_boundary_polyhedron;
      poly_copy<Polyhedron,Exact_polyhedron>(*region_ips[id], *region_eps[id]);
  }

  BOOST_FOREACH(const zone_ep_map::value_type& zone_pair, zone_eps) {
      int id = zone_pair.first;
      Polyhedron* inexact_boundary_polyhedron = new Polyhedron();
      zone_ips[id] = inexact_boundary_polyhedron;
      poly_copy<Polyhedron,Exact_polyhedron>(*zone_ips[id], *zone_eps[id]);
      zone_pips[id] = new Point_inside_polyhedron(*zone_ips[id]);
  }

  if (zone_pips.size())
      std::cout << "Order of zones by priority:";
  BOOST_FOREACH(zone_pip_map::value_type& v, zone_pips) {
      std::cout << " " << v.first;
  }
  if (zone_pips.size())
      std::cout << std::endl;

  Tree *boundary_tree = NULL;
  if (include_boundary_tree) {
      boundary_tree = new Tree();
      //boundary_tree->insert(region_ips[extent]->facets_begin(), region_ips[extent]->facets_end());
      boundary_tree->insert(region_ips[organ]->facets_begin(), region_ips[organ]->facets_end());
      BOOST_FOREACH(std::vector<int>::value_type& id, vessels) {
          boundary_tree->insert(region_ips[id]->facets_begin(), region_ips[id]->facets_end());
      }
  }

  Implicit_zone_function<K> izf(region_ips, centre, bounding_radius, (include_organ ? &organ : NULL), &extent, vessels, needles, zone_pips,
          zone_activity_spheres, default_zone_id);

  FT bbox_radius = 1.01 * pow(CGAL::squared_distance(bbox[0], bbox[7]), .5) / 2;
  if (bbox_radius < 1e-12) {
      std::cerr << "Extent must have non-zero diameter" << std::endl;
      abort();
  }
  else if (bbox_radius > bounding_radius) {
      if (settings.has_bounding_radius())
          std::cout << "(1% padded) Radius of extent bounding box larging than requested bounding radius, using the bounding box radius" << std::endl;
  }

  Mesh_domain_implicit domain(izf, K::Sphere_3(*centre, (FT)1.2*bbox_radius*bbox_radius), (FT)(bbox_radius * 2.e-5));

  Tree *needle_tree = NULL;
  if (include_needle) {
      needle_tree = new Tree();
      BOOST_FOREACH(std::vector<int>::value_type& id, needles) {
          Polyhedron* inexact_needle_polyhedron = new Polyhedron();
          region_ips[id] = inexact_needle_polyhedron;
          poly_copy<Polyhedron,Exact_polyhedron>(*inexact_needle_polyhedron, *region_eps[id]);

          needle_tree->insert(inexact_needle_polyhedron->facets_begin(), inexact_needle_polyhedron->facets_end());
      }
  }
  zone_tree_map zone_trees;
  if (!zone_ips.empty() && !omit_zone_tree) {
      BOOST_FOREACH(const zone_ip_map::value_type& zone_pair, zone_ips) {
          int id = zone_pair.first;
          zone_trees[id] = new Tree();
          //Polyhedron* inexact_zone_polyhedron = new Polyhedron();
          //region_ips[id] = inexact_zone_polyhedron;
          //poly_copy<Polyhedron,Exact_polyhedron>(*inexact_zone_polyhedron, *zone_eps[id]);

          zone_trees[id]->insert(zone_ips[id]->facets_begin(), zone_ips[id]->facets_end());
      }
  }
  Proximity_domain_field_3<Mesh_domain::R,Mesh_domain::Index> pdf(boundary_tree, nearfield, farfield, zone_cls, zone_radius, centre_radius,
          (include_centre ? centre : NULL), (omit_needle_tree ? NULL : needle_tree), zone_trees, granularity, bbox,
          (solid_zone ? &zone_pips : NULL));

  Mesh_criteria new_criteria(cell_size=pdf);

  std::cout << "Tetrahedralizing..." << std::endl;
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, new_criteria, features(domain), no_lloyd(), no_odt(), perturb(10000,1e-3));

  std::cout << " Done" << std::endl;

  if (output_medit) {
      std::cout << "Write to MEDIT" << std::endl;
      std::ofstream file_medit((output_prefix + ".msh").c_str());
      c3t3.output_to_medit(file_medit);
      file_medit.close();
  }

  C3t3::Facets_in_complex_iterator fit, end;
  fit = c3t3.facets_in_complex_begin();
  end = c3t3.facets_in_complex_end();

  /*
  std::cout << "Accelerating distance queries..." << std::endl;
  int k;
  //double dist, dist0, dist1, dist2;

  //if (include_needle)
  //    patches++;
  //Exact_polyhedron* patch_polyhedra[patches+1];
  //patch_polyhedra[1] = &extent; patch_polyhedra[2] = &organ; patch_polyhedra[3] = &pv; patch_polyhedra[4] = &hv; patch_polyhedra[5] = &ha;
  //if (include_needle)
  //    patch_polyhedra[6] = &needle;

  region_tree_map region_trees;
  BOOST_FOREACH(region_ep_map::value_type& region_pair, region_eps) {
      int id = region_pair.first;
      Exact_polyhedron& region = *region_pair.second;

      if (id < 0)
          continue;

      region_trees[id] = new Exact_Tree(region.facets_begin(), region.facets_end());
      region_trees[id]->accelerate_distance_queries();
  }
  */

  //std::string visfile("vis.vts");
  //visualization_save(visfile);

  std::cout << "Assigning boundary indices..." << std::endl;
  Point vi[3];
  Exact_Point v[3];
  K_to_EK conv;
  int l = 0, m = 0;
  std::map<C3t3::Triangulation::Facet,int> boundary_indices;
  region_boundary_count_map region_boundary_count;
  //double length_scale = (bounding_radius < bbox_radius) ? bounding_radius : bbox_radius;
  //std::cout << length_scale << std::endl;
  for( ; fit != end ; ++fit)
  {
      l++;
      if (l%1000 == 0) std::cout << " " << l << "/" << c3t3.number_of_facets_in_complex() << std::endl;

      C3t3::Surface_patch_index spi = c3t3.surface_patch_index(*fit);

      int id = -1;

      /* This facet is an extra-domain boundary!? */
      if (spi.first <= 0 && spi.second <= 0)
          continue;
      /* This facet is on the domain boundary */
      else if (spi.first <= 0)
          id = -spi.first;
      else if (spi.second <= 0)
          id = -spi.second;
      /* This facet is an internal boundary */
      else if (mark_zone_boundaries)
      {
          /* Ensures ordering by sorted priority (this assumes zones and regions are aligned) */
          BOOST_FOREACH(std::vector<int>::value_type& zone_id, zones) {
              if (spi.first == zone_id || spi.second == zone_id)
              {
                  id = zone_id;
                  break;
              }
          }
      }

      if (id < 0)
      {
          boundary_indices[(*fit)] = 0;
          continue;
      }

      m++;

      boundary_indices[(*fit)] = id;

      if (!region_boundary_count.count(id))
          region_boundary_count[id] = 0;
      region_boundary_count[id]++;
  }
  std::cout << m << " / " << l << " surface facets tagged" << std::endl;
  BOOST_FOREACH(region_boundary_count_map::value_type& region_pair, region_boundary_count) {
      std::cout << " of which " << region_pair.second << " in region " << region_pair.first << std::endl;
  }


  //C3t3::Cells_in_complex_iterator cit, cend;
  //cit = c3t3.cells_in_complex_begin();
  //cend = c3t3.cells_in_complex_end();

  //std::cout << "Assigning internal indices..." << std::endl;
  //l = 0; m = 0;
  //std::map<C3t3::Triangulation::Cell_handle,int> domain_indices;
  //for( ; cit != cend ; ++cit)
  //{
  //    l++;
  //    if (l%1000 == 0) std::cout << " " << l << "/" << c3t3.number_of_cells_in_complex() << std::endl;
  //    int n = 0;

  //    BOOST_FOREACH(zone_pip_map::value_type& zone_pair, zone_pips) {
  //        int id = zone_pair.first;
  //        Point_inside_polyhedron& zone = *(zone_pair.second);

  //        if (zone(CGAL::circumcenter(c3t3.triangulation().tetrahedron(cit)))) {
  //            domain_indices[cit] = id;
  //            break;
  //        }
  //        else {
  //            domain_indices[cit] = 0;
  //        }
  //   }
  //}
  //std::cout << m << " / " << l << " volume tetrahedra tagged" << std::endl;

  if (output_gmsh) {
      std::cout << "Write to GMSH" << std::endl;
      CGAL::write_c3t3_to_gmsh_file<C3t3,Polyhedron>(c3t3, boundary_indices, output_prefix + ".msh", (default_zone_id == 0), number_from_zero);
  }
  if (output_vtk) {
      //Removed until licencing clarified.
      //std::cout << "Write to VTK" << std::endl;
      //CGAL::write_c3t3_to_vtk_xml_file<C3t3>(c3t3, output_prefix + ".vtu");
  }

  std::cout << "Complete" << std::endl;
  return 0;
}
