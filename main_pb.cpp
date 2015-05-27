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

#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>


using namespace mesherCGAL;

/* Convenience standalone to allow compilation and execution without any non-header boost libraries, i.e. program_options */
int main(int argc, char* argv[])
{
  if (argc < 2 || argc > 4) {
      std::cerr << "There should be exactly one argument representing the settings file" << std::endl;
      exit(-2);
  }

  CGALSettings settings;

  std::ifstream from_file(argv[1]);
  //std::cout << argv[2] << std::endl;
  if (argc > 2 && strcmp(argv[2], "-text") == 0) {
	  //std::cout << argv[1] << std::endl;
	  //while (!from_file.eof())
	  //{
		 // std::string fromeee;
		 // getline(from_file, fromeee);
		 // std::cout <<"hey"<< fromeee <<std::endl;
		 // 
	  //}
      google::protobuf::io::IstreamInputStream is(&from_file);
      if (!google::protobuf::TextFormat::Parse(&is, &settings)) {
          std::cerr << "Unable to parse protobuf from the settings file (text)" << std::endl;
          exit(-3);
      }
  } else {
      if (!settings.ParseFromIstream(&from_file)) {
          std::cerr << "Unable to parse protobuf from the settings file (binary)" << std::endl;
          exit(-3);
      }
  }

  from_file.close();

  if (argc > 2 && strcmp(argv[2], "-print") == 0) {
      std::string settings_string;
      google::protobuf::TextFormat::PrintToString(settings, &settings_string);

      std::cout << "LOADED SETTINGS FROM PROTOBUF" << std::endl
          << "+++++++++++++++++++++++++++++++" << std::endl
          << settings_string
          << "+++++++++++++++++++++++++++++++" << std::endl;

      return 0;
  }

  if (argc > 2 && strcmp(argv[2], "-dummy") == 0)
      return 0;

  mesherCGAL::run(settings);
}


