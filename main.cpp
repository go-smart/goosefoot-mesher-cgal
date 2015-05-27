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

#include "mesher_cgal.h"

#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
namespace cliopt = boost::program_options;

using namespace mesherCGAL;

void split_region(const std::string& input_string, std::string& file_string, int& index, float& characteristic_length, float& priority) {
  std::vector<std::string> pair;
  boost::split(pair, input_string, boost::is_any_of(":"));

  index = boost::lexical_cast<int>(pair[1]);
  file_string = pair[0];

  if (pair.size() > 2)
      characteristic_length = boost::lexical_cast<float>(pair[2]);
  if (pair.size() > 3)
      priority = boost::lexical_cast<float>(pair[3]);
  else
      priority = 0;
}

int main(int argc, char* argv[])
{
  bool tetrahedralize_only = false, output_medit = false, output_gmsh = true, output_vtk = true, verbose = false, suppress_nef_loading = false,
       mark_zone_boundaries = false, solid_zone = false;
  std::string organ_file("cgal/organ.off"),
              ha_file("cgal/ha.off"),
              pv_file("cgal/pv.off"),
              extent_file("cgal/extent.off"),
              needle_file("cgal/needle.off"),
              combined_file("cgal/_combined.off"),
              boundary_file("cgal/_boundary.off"),
              structures_file("cgal/_structure.off"),
              needles_string,
              vessels_string,
              zones_string;
  float centre_radius;
  float zone_radius;
  float bounding_radius;
  int tissue_id;
  bool dump_settings;
  std::string settings_binary_file;

  cliopt::options_description options_description("Allowed options");
  options_description.add_options()
    ("help,h,?", "produce help message")
    ("dump_settings,D", cliopt::value(&dump_settings)->zero_tokens(), "dump settings Protocol Buffer (as ASCII)")
    ("dump_settings_binary", cliopt::value<std::string>(&settings_binary_file), "dump settings Protocol Buffer (as ASCII)")

    ("centre,c", cliopt::value< std::vector<float> >()->multitoken(), "set conceptual centre of mesh")
    ("centre_radius,r", cliopt::value<float>(&centre_radius)->default_value(1.), "set inner radius around centre for high density cells")
    ("zone_radius", cliopt::value<float>(&zone_radius)->default_value(0.), "set radius around zone for high density cells")
    ("tissueid,t", cliopt::value<int>(&tissue_id)->default_value(1), "default ID of volume elements not otherwise in a zone")
    ("bounding_radius,R", cliopt::value<float>(&bounding_radius)->default_value(50.), "set radius around extent (NB: overridden if extent dimensions exceed this)")
    ("tetrahedralize-only,t", "do not combine surfaces but use previously produced combined surfaces for tetrahedralization")
    ("suppress-nef-loading,s", cliopt::value(&suppress_nef_loading)->zero_tokens(), "suppress loading of cached NEF polyhedra if newer than input")
    ("mark_zone_boundaries,M", cliopt::value(&mark_zone_boundaries)->zero_tokens(), "include internal boundaries in the output mesh")
    ("solid_zone", cliopt::value(&solid_zone)->zero_tokens(), "make inside of zone uniform minimum zone_cl characteristic length")
    ("boundary_tree,b", "include boundary tree in distance calculations")
    ("omit_needle_tree,d", "do not include needle tree in distance calculations")
    ("omit_zone_tree,q", "do not include zone tree in distance calculations")
    ("dense_centre", "focus a region of high density cells on the centre")
    ("number_from_zero", "start the point numbering from zero")
    ("matching_tolerance", "how tightly the surfaces are expected to match their meshed surfaces")

    ("nearfield,n", cliopt::value<float>()->default_value(0.3), "near-field characteristic length")
    ("farfield,f", cliopt::value<float>()->default_value(3.), "far-field characteristic length")
    ("zonefield,z", cliopt::value<float>()->default_value(0.3), "characteristic length near zone boundaries")
    ("granularity,g", cliopt::value<float>()->default_value(1.), "granularity (dimension) of caching blocks")

    ("organ,L", cliopt::value<std::string>(&organ_file), "location of organ surface")
    ("hepatic_artery,H", cliopt::value<std::string>(&ha_file), "location of hepatic artery surface")
    ("vessels,V", cliopt::value< std::vector<std::string> >()->multitoken(), "location of vessel surfaces")
    ("zones,Z", cliopt::value< std::vector<std::string> >()->multitoken(), "location of internal zones")
    ("portal_vein,P", cliopt::value<std::string>(&pv_file), "location of hepatic portal vein surface")
    ("extent,S", cliopt::value<std::string>(&extent_file), "location of extent surface")
    ("needles,N", cliopt::value< std::vector<std::string> >()->multitoken(), "location of needle surfaces")
    ("combined_surface,C", cliopt::value<std::string>(&combined_file), "input/output location of combined surface")
    ("boundary_surface,B", cliopt::value<std::string>(&boundary_file), "input/output location of outer (boundary) surface")
    ("structure_surface,I", cliopt::value<std::string>(&structures_file), "input/output location of internal structure surface")

    ("output_prefix,p", cliopt::value<std::string>(), "prefix of the output mesh")
    ("output_vtk,k", "output VTK")
    ("output_gmsh,m", "output MSH for GMSH")
    ("output_medit", "output Medit")

    ("number-from-zero", "number elements from zero in MSH output")

    ("verbose,v", "verbose output")
    ("version", "help and version")
  ;

  cliopt::variables_map vm;
  cliopt::store(cliopt::parse_command_line(argc, argv, options_description, cliopt::command_line_style::unix_style ^ cliopt::command_line_style::allow_short), vm);
  cliopt::notify(vm);

  CGALSettings settings;

  settings.set_output_medit(!vm["output_medit"].empty());
  settings.set_number_from_zero(!vm["number_from_zero"].empty());
  settings.set_verbose(!vm["verbose"].empty());
  settings.set_output_gmsh(!vm["output_gmsh"].empty());
  settings.set_output_vtk(!vm["output_vtk"].empty());
  settings.set_dense_centre(!vm["dense_centre"].empty());
  settings.set_omit_needle_tree(!vm["omit_needle_tree"].empty());
  settings.set_boundary_tree(!vm["boundary_tree"].empty());
  settings.set_omit_zone_tree(!vm["omit_zone_tree"].empty());
  settings.set_mark_zone_boundaries(!vm["mark_zone_boundaries"].empty());
  settings.set_solid_zone(!vm["solid_zone"].empty());
  settings.set_dump_settings(dump_settings);

  if (!vm["nearfield"].empty())
      settings.set_near_field(vm["nearfield"].as<float>());
  if (!vm["farfield"].empty())
      settings.set_far_field(vm["farfield"].as<float>());
  if (!vm["zonefield"].empty())
      settings.set_zone_field(vm["zonefield"].as<float>());
  if (!vm["granularity"].empty())
      settings.set_granularity(vm["granularity"].as<float>());
  if (!vm["output_prefix"].empty())
      settings.set_output_prefix(vm["output_prefix"].as<std::string>());

  if (vm.count("help")) {
    std::cout << "CGAL Mesher - 0.1" << std::endl;
    std::cout << options_description << std::endl;
    return 1;
  }

  settings.set_centre_radius(centre_radius);
  settings.set_tissue_id(tissue_id);
  settings.set_bounding_radius(bounding_radius);
  settings.set_tetrahedralize_only(!vm["tetrahedralize-only"].empty());

  std::vector<float> coords;
  if (vm.count("centre")) {
      coords = vm["centre"].as< std::vector<float> >();
      if (coords.size() != 3) {
          std::cerr << "Need an \"x,y,z\" argument to define centre point - need 3 coords not " << coords.size() << std::endl;
          abort();
      }

      settings.add_centre(coords[0]);
      settings.add_centre(coords[1]);
      settings.add_centre(coords[2]);
  }

  std::string file_string;
  int index;
  float characteristic_length, priority;

  bool include_organ = (!vm["organ"].empty());

  if (include_organ) {
       split_region(organ_file, file_string, index, characteristic_length, priority);
       settings.set_organ_file(file_string);
       settings.set_organ_index(index);
  }

  if (!vm["extent"].empty()) {
      split_region(extent_file, file_string, index, characteristic_length, priority);
      settings.set_extent_file(file_string);
      settings.set_extent_index(index);
  }

  if (vm.count("needles")) {
      BOOST_FOREACH(const std::string& needle_string, vm["needles"].as< std::vector< std::string > >()) {
          split_region(needle_string, file_string, index, characteristic_length, priority);
          settings.add_needles(file_string);
          settings.add_needles_index(index);
      }
  }

  if (vm.count("vessels")) {
      BOOST_FOREACH(const std::string& vessel_string, vm["vessels"].as< std::vector< std::string > >()) {
          split_region(vessel_string, file_string, index, characteristic_length, priority);
          settings.add_vessels(file_string);
          settings.add_vessels_index(index);
      }
  }

  if (vm.count("zones")) {
      BOOST_FOREACH(const std::string& zone_string, vm["zones"].as< std::vector< std::string > >()) {
          float characteristic_length = -1, priority = 0;
          split_region(zone_string, file_string, index, characteristic_length, priority);
          CGALSettings::Zone* zone = settings.add_zone();
          zone->set_file(file_string);
          zone->set_index(index);
          zone->set_priority(priority);
          if (characteristic_length > 0)
              zone->set_characteristic_length(characteristic_length);
      }
  }

  if (!vm["dump_settings_binary"].empty()) {
    std::ofstream temp_settings_fd(settings_binary_file.c_str(), std::ios_base::out | std::ios_base::binary);

    if (!settings.SerializeToOstream(&temp_settings_fd)) {
        std::cerr << "Could not output settings file" << std::endl;
    }

    temp_settings_fd.close();
  }

  //if (!vm["combined_surface"].empty()) {
  //    combined = parse_region(region_files, combined_file);
  //    std::cout << " - Combined : " << combined << std::endl;
  //}

  //if (!vm["boundary_surface"].empty()) {
  //    boundary_id = parse_region(region_files, boundary_file);
  //    std::cout << " - Boundary : " << boundary_id << std::endl;
  //}

  //if (!vm["structures_surfaces"].empty()) {
  //    structures = parse_region(region_files, structures_file);
  //    std::cout << " - Structures : " << structures << std::endl;
  //}

  mesherCGAL::run(settings);
}


