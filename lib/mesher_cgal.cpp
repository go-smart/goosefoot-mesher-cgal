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

#include "mesher_cgal_app.h"

// [B]
//
// Tree
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
// ]A]
#include <CGAL/Side_of_triangle_mesh.h>
// [C]

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

#include <vtkXMLUnstructuredGridReader.h>

#include <map>
#include <vector>
#include <utility>

#include "parse_region.h"
#include "write_c3t3_to_gmsh_file.h"
#include "zone_priority_sorting.h"
#include "copy_polyhedron.h"

// [B]
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
// End [B]

// [A]
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

typedef CGAL::AABB_face_graph_triangle_primitive<Exact_polyhedron> Exact_Primitive;
typedef CGAL::AABB_traits<Exact_Kernel, Exact_Primitive> Exact_Traits;
typedef CGAL::AABB_tree<Exact_Traits> Exact_Tree;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

typedef std::map< int, int > region_boundary_count_map;
typedef std::map< int, std::string > region_string_map;
typedef std::map< int, Exact_Tree* > region_tree_map;

typedef std::map< int, std::string > zone_string_map;
typedef std::map< int, Exact_polyhedron* > zone_ep_map;
// End [A]
// [C]
typedef CGAL::Side_of_triangle_mesh<Polyhedron,K> Side_of_triangle_mesh;
// End[C]

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Exact_Mesh_criteria;

bool mesherCGAL::Zone::set_activity_sphere(float x, float y, float z, float r, int i) {
    activity_sphere = new struct activity_sphere;
    activity_sphere->centre = new Point(x, y, z);
    activity_sphere->r = r;
    activity_sphere->i = i;
    return true;
}

void mesherCGAL::Zone::print_activity_sphere() {
    if (activity_sphere)
      std::cout << "[" << activity_sphere->i <<
          ":(" << activity_sphere->centre->x() <<
          "," << activity_sphere->centre->y() <<
          "," << activity_sphere->centre->z() <<
          "):" << activity_sphere->r << "]";
}

bool mesherCGAL::PolyhedralZone::load(std::string filename) {
    std::shared_ptr<Exact_polyhedron> ep = std::make_shared<Exact_polyhedron>();
    PolyhedronUtils::readSurfaceFile(filename, *ep);
    create_polyhedra(ep);

    return true;
}

bool mesherCGAL::PolyhedralZone::create_polyhedra(std::shared_ptr<Exact_polyhedron> ep) {
    _ep = ep;

    _ip = std::make_shared<Polyhedron>();
    poly_copy<Polyhedron,Exact_polyhedron>(_ip, _ep);
    _pip = std::make_shared<Side_of_triangle_mesh>(*_ip);

    return true;
}

bool mesherCGAL::UnstructuredGridZone::load(std::string filename) {
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

    reader->SetFileName(filename.c_str());
    reader->Update();
    _grid = reader->GetOutput();

    return true;
}

bool mesherCGAL::UnstructuredGridZone::contains_all(const K::Point_3& p) const {
    double x[3] = {p.x(), p.y(), p.z()}, pcoords[3], weights[4];
    int subId;
    vtkIdType cid = _grid->FindCell(x, NULL, 0, 1e-10, subId, pcoords, weights);
    return cid < 0 ? 0 : _id;
}

int mesherCGAL::MesherCGAL::init() {
  std::string _settings_string;
  google::protobuf::TextFormat::PrintToString(_settings, &_settings_string);
  if (_settings.dump_settings()) {
      std::cout << "LOADED SETTINGS FROM PROTOBUF" << std::endl
          << "+++++++++++++++++++++++++++++++" << std::endl
          << _settings_string
          << "+++++++++++++++++++++++++++++++" << std::endl << std::endl;
  }
  std::vector<float> coords;

  if (_settings.centre_size() > 0) {
      if (_settings.centre_size() != 3) {
          std::cerr << "ERROR: Need an \"x,y,z\" argument to define centre point - need 3 coords not " << coords.size() << std::endl;
          exit(4);
      }

      _centre = new Point(_settings.centre(0), _settings.centre(1), _settings.centre(2));

      std::cout << "Centre is set at (" << _centre->x() << ", " << _centre->y() << ", " << _centre->z() << ")" << std::endl;
  }
  if (!_settings.has_granularity())
      _settings.set_granularity(1.);
  if (!_settings.has_bounding_radius())
      _settings.set_bounding_radius(50.);

  return 0;
}

int mesherCGAL::MesherCGAL::setup_regions() {
  region_string_map region_files;

  std::cout << "Finding regions : " << std::endl;
  if (_settings.has_organ_file()) {
      region_files[_settings.organ_index()] = _settings.organ_file();
      std::cout << " - Organ : " << _settings.organ_index() << std::endl;
  }
  std::cout << " - Extent : " << _settings.extent_index() << std::endl;
  if (_settings.needles_size() > 0) {
      std::cout << " - Needles : " << std::endl << "   \\ ";
      for (int i = 0 ; i < _settings.needles_size() ; i++) {
          int id = _settings.needles_index(i);
          region_files[id] = _settings.needles(i);
          needles.push_back(id);
          std::cout << id << " ";
      }
      std::cout << std::endl;
  }
  if (_settings.vessels_size()) {
      std::cout << " - Vessels : " << std::endl << "   \\ ";
      for (int i = 0 ; i < _settings.vessels_size() ; i++) {
          int id = _settings.vessels_index(i);
          region_files[id] = _settings.vessels(i);
          vessels.push_back(id);
          std::cout << id << " ";
      }
      std::cout << std::endl;
  }
  if (_settings.zones_size() || _settings.zone_size()) {
      std::cout << " - Zones : " << std::endl << "   \\ ";
      for (int i = 0 ; i < _settings.zone_size() ; i++) {
          int id = _settings.zone(i).index();
          float cl, priority;
          std::cout << id;

          if (_settings.zone(i).has_characteristic_length()) {
              cl = _settings.zone(i).characteristic_length();
              std::cout << "(" << cl << ")";
          }
          else {
              cl = _settings.zone_field();
          }

          if (_settings.zone(i).has_priority()) {
              priority = _settings.zone(i).priority();
              if (fabs(priority) > 1e-10)
                  std::cout << "<" << priority << ">";
          }
          else {
              priority = 0.0;
          }

          std::string filename = _settings.zone(i).file();

          if (filename.substr(std::max(4, (int)filename.length()) - 4) == std::string(".vtu")) {
              _zones.push_back(std::unique_ptr<UnstructuredGridZone>(new UnstructuredGridZone(id, cl, priority)));
              _zones.back()->load(filename);
          }
          else {
              /* In lieu of make_unique in C++14 */
              _zones.push_back(std::unique_ptr<PolyhedralZone>(new PolyhedralZone(id, cl, priority)));
              Zone& zone = *_zones.back();

              zone.load(filename);

              std::shared_ptr<Exact_polyhedron> ep;
              if ((ep = zone.exact_polyhedron())) {
                  if (_settings.has_mark_zone_boundaries() && _settings.mark_zone_boundaries())
                      _region_eps[id] = ep;

                  std::cout << id << " (" << ep->size_of_facets() << " facets) " << std::flush;
              }
          }
          Zone& zone = *_zones.back();

          if (_settings.zone(i).has_activity() && _settings.zone(i).has_inactivity_index()) {
              double
                  x = _settings.zone(i).activity().x(),
                  y = _settings.zone(i).activity().y(),
                  z = _settings.zone(i).activity().z();
              if (_centre != NULL)
              {
                  x += _centre->x();
                  y += _centre->y();
                  z += _centre->z();
              }
              //int inactivity_index = _settings.zone(i).inactivity_index();
              zone.set_activity_sphere(x, y, z, _settings.zone(i).activity().r(), i);
              //_zones.push_back(PolyhedralZone(inactivity_index, cl, priority));
              zone.print_activity_sphere();
          }
          std::cout << " ";
      }
      std::cout << std::endl;
  }
  std::sort(_zones.begin(), _zones.end(), ZonePrioritySorting());

  Polyhedron inexact_structures_polyhedron, inexact_combined_polyhedron, inexact_boundary_polyhedron;
  Exact_polyhedron exact_structures_polyhedron;
  Exact_polyhedron exact_boundary_polyhedron;
  Exact_polyhedron exact_combined_polyhedron;
  std::cout << "Loading (boundary) regions : " << std::endl << " - ";
  BOOST_FOREACH(const region_string_map::value_type& region_pair, region_files) {
      std::shared_ptr<Exact_polyhedron> ep = std::make_shared<Exact_polyhedron>();
      PolyhedronUtils::readSurfaceFile(region_pair.second, *ep);
      _region_eps[region_pair.first] = ep;
      std::cout << region_pair.first << " (" << ep->size_of_facets() << " facets) " << std::flush;
  }
  std::cout << std::endl;
  std::cout << std::endl;
  if (_settings.has_organ_file()) {
      std::shared_ptr<Polyhedron> inexact_organ_polyhedron = std::make_shared<Polyhedron>();
      _region_ips[_settings.organ_index()] = inexact_organ_polyhedron;
      poly_copy<Polyhedron,Exact_polyhedron>(inexact_organ_polyhedron, _region_eps[_settings.organ_index()]);
  }

  BOOST_FOREACH(std::vector<int>::value_type& id, vessels) {
      std::shared_ptr<Polyhedron> inexact_boundary_polyhedron = std::make_shared<Polyhedron>();
      _region_ips[id] = inexact_boundary_polyhedron;
      poly_copy<Polyhedron,Exact_polyhedron>(_region_ips[id], _region_eps[id]);
  }
  BOOST_FOREACH(std::vector<int>::value_type& id, needles) {
      std::shared_ptr<Polyhedron> inexact_boundary_polyhedron = std::make_shared<Polyhedron>();
      _region_ips[id] = inexact_boundary_polyhedron;
      poly_copy<Polyhedron,Exact_polyhedron>(_region_ips[id], _region_eps[id]);
  }

  if (_zones.size())
      std::cout << "Order of zones by priority:";
  for (auto&& v : _zones)
      std::cout << " " << v->get_id();
  if (_zones.size())
      std::cout << std::endl;

  if (_settings.boundary_tree()) {
      _boundary_tree = new Tree();
      _boundary_tree->insert(_region_ips[_settings.organ_index()]->facets_begin(), _region_ips[_settings.organ_index()]->facets_end(), *_region_ips[_settings.organ_index()]);
      BOOST_FOREACH(std::vector<int>::value_type& id, vessels) {
          _boundary_tree->insert(_region_ips[id]->facets_begin(), _region_ips[id]->facets_end(), *_region_ips[id]);
      }
  }

  return 0;
}

int mesherCGAL::MesherCGAL::calculate_bbox() {
  if (_settings.has_organ_file()) {
      _bbox_p = new Iso_cuboid(CGAL::bounding_box(_region_ips[_settings.organ_index()]->points_begin(), _region_ips[_settings.organ_index()]->points_end()));
      if (_centre == NULL)
          _centre = new Point(((*_bbox_p)[0].x() + (*_bbox_p)[7].x()) / 2,
                             ((*_bbox_p)[0].y() + (*_bbox_p)[7].y()) / 2,
                             ((*_bbox_p)[0].z() + (*_bbox_p)[7].z()) / 2);
      Side_of_triangle_mesh organ_pip(*_region_ips[_settings.organ_index()]);
      if (organ_pip(*_centre) == CGAL::ON_UNBOUNDED_SIDE)
      {
          std::cerr << "ERROR: Centre lies outside organ" << std::endl;
          return 1;
      }
  } else if (_centre != NULL) {
      _bbox_p = new Iso_cuboid(
          _centre->x() - _settings.bounding_radius(),
          _centre->y() - _settings.bounding_radius(),
          _centre->z() - _settings.bounding_radius(),
          _centre->x() + _settings.bounding_radius(),
          _centre->y() + _settings.bounding_radius(),
          _centre->z() + _settings.bounding_radius()
      );
  } else {
      std::cerr << "ERROR: No centre defined and no organ to guess from" << std::endl;
      return 2;
  }
  const Iso_cuboid& bbox(*_bbox_p);
  _bbox_radius = 1.01 * pow(CGAL::squared_distance(bbox[0], bbox[7]), .5) / 2;
  if (_bbox_radius < 1e-12) {
      std::cerr << "ERROR: Extent must have non-zero diameter" << std::endl;
      return 3;
  }
  else if (_bbox_radius > _settings.bounding_radius()) {
      std::cout << "(1% padded) Radius of extent bounding box larging than requested bounding radius, using the bounding box radius" << std::endl;
  }
   return 0;
}

int mesherCGAL::MesherCGAL::setup_domain_field() {
  Implicit_zone_function<K> izf(_region_ips, _centre, _settings.bounding_radius(), _settings.has_organ_file(), _settings.organ_index(), _settings.has_extent_index(), _settings.extent_index(), vessels, needles, _zones,
          _settings.tissue_id());
  _domain = new Mesh_domain_implicit(izf, K::Sphere_3(*_centre, (FT)1.2*_bbox_radius*_bbox_radius), (FT)(_bbox_radius * 2.e-5));
  Tree *needle_tree = NULL;
  if (_settings.needles_size() > 0) {
      needle_tree = new Tree();
      BOOST_FOREACH(std::vector<int>::value_type& id, needles) {
          std::shared_ptr<Polyhedron> inexact_needle_polyhedron = std::make_shared<Polyhedron>();
          _region_ips[id] = inexact_needle_polyhedron;
          poly_copy<Polyhedron,Exact_polyhedron>(inexact_needle_polyhedron, _region_eps[id]);

          needle_tree->insert(inexact_needle_polyhedron->facets_begin(), inexact_needle_polyhedron->facets_end(), *inexact_needle_polyhedron);
      }
  }
  if (!_settings.omit_zone_tree())
      for (auto&& z : _zones)
          z->add_tree();

  _pdf = new Proximity_domain_field_3<Mesh_domain::R,Mesh_domain::Index>(_boundary_tree, _settings.near_field(), _settings.far_field(), _zones, _settings.zone_radius(), _settings.centre_radius(),
          (_settings.dense_centre() ? _centre : NULL), (_settings.omit_needle_tree() ? NULL : needle_tree), _settings.granularity(), *_bbox_p,
          _settings.has_solid_zone() && _settings.solid_zone());

  return 0;
}

int mesherCGAL::MesherCGAL::mesh() {
  setup_domain_field();
  Mesh_criteria new_criteria(cell_size=*_pdf);

  std::cout << "Tetrahedralizing..." << std::endl;
  _c3t3 = CGAL::make_mesh_3<C3t3>(*_domain, new_criteria, features(*_domain), no_lloyd(), no_odt(), perturb(10000,1e-3));

  std::cout << " Done" << std::endl;

  std::map< int, int > zone_count;
  for (C3t3::Cells_in_complex_iterator it = _c3t3.cells_in_complex_begin(); it != _c3t3.cells_in_complex_end(); ++it)
  {
      int subId = _c3t3.subdomain_index(it);
      if (zone_count.find(subId) == zone_count.end())
          zone_count[subId] = 0;
      zone_count[subId]++;
  }
  for (auto&& pair : zone_count) {
      int subId = pair.first;
      if (_settings.tissue_id() == 0)
          subId--;
      std::cout << pair.second << " cells with id " << subId << std::endl;
  }

  return 0;
}

int mesherCGAL::MesherCGAL::label_boundaries() {
  C3t3::Facets_in_complex_iterator fit, end;
  fit = _c3t3.facets_in_complex_begin();
  end = _c3t3.facets_in_complex_end();
  std::cout << "Assigning boundary indices..." << std::endl;
  Exact_Point v[3];
  int l = 0, m = 0;
  region_boundary_count_map region_boundary_count;
  for( ; fit != end ; ++fit)
  {
      l++;
      if (l%1000 == 0) std::cout << " " << l << "/" << _c3t3.number_of_facets_in_complex() << std::endl;

      C3t3::Surface_patch_index spi = _c3t3.surface_patch_index(*fit);

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
      else if (_settings.has_mark_zone_boundaries() && _settings.mark_zone_boundaries())
      {
          /* Ensures ordering by sorted priority (this assumes zones and regions are aligned) */
          for (auto&& z : _zones)
              if (spi.first == z->get_id() || spi.second == z->get_id())
              {
                  id = z->get_id();
                  break;
              }
      }

      if (id < 0)
      {
          _boundary_indices[(*fit)] = 0;
          continue;
      }

      m++;

      _boundary_indices[(*fit)] = id;

      if (!region_boundary_count.count(id))
          region_boundary_count[id] = 0;
      region_boundary_count[id]++;
  }
  std::cout << m << " / " << l << " surface facets tagged" << std::endl;
  BOOST_FOREACH(region_boundary_count_map::value_type& region_pair, region_boundary_count) {
      std::cout << " of which " << region_pair.second << " in region " << region_pair.first << std::endl;
  }

  return 0;
}

int mesherCGAL::MesherCGAL::output() {
  std::cout << "Write to GMSH" << std::endl;
  CGAL::write_c3t3_to_gmsh_file<C3t3,Polyhedron>(_c3t3, _boundary_indices, std::string(_settings.output_prefix()) + ".msh", (_settings.tissue_id() == 0), _settings.number_from_zero());

  std::cout << "Complete" << std::endl;

  return 0;
}

int mesherCGAL::MesherCGAL::run() {
    int err;

    if ((err = setup_regions())) {
        std::cerr << "Could not set up regions" << std::endl;
        return 1000 + err;
    }
    if ((err = calculate_bbox())) {
        std::cerr << "Could not calculate BBox" << std::endl;
        return 2000 + err;
    }

    if ((err = mesh())) {
        std::cerr << "Could not mesh" << std::endl;
        return 4000 + err;
    }

    if ((err = label_boundaries())) {
        std::cerr << "Could not label boundaries" << std::endl;
        return 5000 + err;
    }

    if ((err = output())) {
        std::cerr << "Could not output" << std::endl;
        return 6000 + err;
    }

    return 0;
}
