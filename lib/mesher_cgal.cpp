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

#include "copy_polyhedron.h"
#include "mesher_cgal.h"

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

typedef K::Triangle_3 Triangle;
typedef K::Plane_3 Plane;
typedef K::Vector_3 Vector;
typedef K::Segment_3 Segment;
typedef K::Iso_cuboid_3 Iso_cuboid;

// Tree
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef CGAL::AABB_face_graph_triangle_primitive<Exact_polyhedron> Exact_Primitive;
typedef CGAL::AABB_traits<Exact_Kernel, Exact_Primitive> Exact_Traits;
typedef CGAL::AABB_tree<Exact_Traits> Exact_Tree;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

typedef std::map< int, struct activity_sphere* > zone_activity_sphere_map;
typedef std::map< int, int > region_boundary_count_map;
typedef std::map< int, std::string > region_string_map;
typedef std::map< int, Exact_polyhedron* > region_ep_map;
typedef std::map< int, Polyhedron* > region_ip_map;
typedef std::map< int, Exact_Tree* > region_tree_map;

typedef std::map< int, std::string > zone_string_map;
typedef std::map< int, Polyhedron* > zone_ip_map;
typedef std::map< int, Exact_polyhedron* > zone_ep_map;

#include "zone_priority_sorting.h"

#include "Signed_mesh_domain_3.h"
#include "implicit_zone_function.h"

#include <CGAL/Point_inside_polyhedron_3.h>
typedef CGAL::Point_inside_polyhedron_3<Polyhedron,K> Point_inside_polyhedron;

typedef CGAL::Signed_mesh_domain_3< Implicit_zone_function<K>, K > Mesh_domain_implicit;

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

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

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

#include "mesher_cgal.int.h"

#include <google/protobuf/text_format.h>

#include <CGAL/Simple_cartesian.h>

#include <map>
#include <vector>
#include <utility>

#include "parse_region.h"
#include "write_c3t3_to_gmsh_file.h"

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Exact_Mesh_criteria;

int mesherCGAL::run(CGALSettings& settings) {

  std::string output_prefix("gsmesh.out");
  bool include_boundary_tree = false, omit_needle_tree = false, omit_zone_tree = false, include_needle = false,
       include_centre = false, solid_zone = false, mark_zone_boundaries = false;
  bool number_from_zero = false;
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
          std::cerr << "ERROR: Need an \"x,y,z\" argument to define centre point - need 3 coords not " << coords.size() << std::endl;
          exit(4);
      }

      centre = new Point(settings.centre(0), settings.centre(1), settings.centre(2));

      std::cout << "Centre is set at (" << centre->x() << ", " << centre->y() << ", " << centre->z() << ")" << std::endl;
  }

  number_from_zero = settings.number_from_zero();
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
              std::cout << "[" << zone_activity_spheres[id]->i <<
                  ":(" << zone_activity_spheres[id]->centre->x() <<
                  "," << zone_activity_spheres[id]->centre->y() <<
                  "," << zone_activity_spheres[id]->centre->z() <<
                  "):" << zone_activity_spheres[id]->r << "]";
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

      Point_inside_polyhedron organ_pip(*region_ips[organ]);
      if (!organ_pip(*centre))
      {
          std::cerr << "ERROR: Centre lies outside organ" << std::endl;
          exit(2);
      }
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
      std::cerr << "ERROR: No centre defined and no organ to guess from" << std::endl;
      exit(3);
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
      boundary_tree->insert(region_ips[organ]->facets_begin(), region_ips[organ]->facets_end(), *region_ips[organ]);
      BOOST_FOREACH(std::vector<int>::value_type& id, vessels) {
          boundary_tree->insert(region_ips[id]->facets_begin(), region_ips[id]->facets_end(), *region_ips[id]);
      }
  }

  Implicit_zone_function<K> izf(region_ips, centre, bounding_radius, (include_organ ? &organ : NULL), &extent, vessels, needles, zone_pips,
          zone_activity_spheres, default_zone_id);

  FT bbox_radius = 1.01 * pow(CGAL::squared_distance(bbox[0], bbox[7]), .5) / 2;
  if (bbox_radius < 1e-12) {
      std::cerr << "ERROR: Extent must have non-zero diameter" << std::endl;
      exit(3);
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

          needle_tree->insert(inexact_needle_polyhedron->facets_begin(), inexact_needle_polyhedron->facets_end(), *inexact_needle_polyhedron);
      }
  }
  zone_tree_map zone_trees;
  if (!zone_ips.empty() && !omit_zone_tree) {
      BOOST_FOREACH(const zone_ip_map::value_type& zone_pair, zone_ips) {
          int id = zone_pair.first;
          zone_trees[id] = new Tree();
          zone_trees[id]->insert(zone_ips[id]->facets_begin(), zone_ips[id]->facets_end(), *zone_ips[id]);
      }
  }
  Proximity_domain_field_3<Mesh_domain::R,Mesh_domain::Index> pdf(boundary_tree, nearfield, farfield, zone_cls, zone_radius, centre_radius,
          (include_centre ? centre : NULL), (omit_needle_tree ? NULL : needle_tree), zone_trees, granularity, bbox,
          (solid_zone ? &zone_pips : NULL));

  Mesh_criteria new_criteria(cell_size=pdf);

  std::cout << "Tetrahedralizing..." << std::endl;
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, new_criteria, features(domain), no_lloyd(), no_odt(), perturb(10000,1e-3));

  std::cout << " Done" << std::endl;

  C3t3::Facets_in_complex_iterator fit, end;
  fit = c3t3.facets_in_complex_begin();
  end = c3t3.facets_in_complex_end();

  std::cout << "Assigning boundary indices..." << std::endl;
  Exact_Point v[3];
  int l = 0, m = 0;
  std::map<C3t3::Triangulation::Facet,int> boundary_indices;
  region_boundary_count_map region_boundary_count;

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
          BOOST_FOREACH(const zone_pip_map::value_type& zone_pair, zone_pips) {
              if (spi.first == zone_pair.first || spi.second == zone_pair.first)
              {
                  id = zone_pair.first;
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


  std::cout << "Write to GMSH" << std::endl;
  CGAL::write_c3t3_to_gmsh_file<C3t3,Polyhedron>(c3t3, boundary_indices, output_prefix + ".msh", (default_zone_id == 0), number_from_zero);

  std::cout << "Complete" << std::endl;

  return 0;
}
