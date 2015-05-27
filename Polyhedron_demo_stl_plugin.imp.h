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

/*
 * This file is adapted under the very permissive CGAL-4.3 examples licence
 * Note that there were therefore unnamed upstream authors. The license
 * is reproduced as follows:
 *     Copyright (c) 1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007
 *     Utrecht University (The Netherlands),
 *     ETH Zurich (Switzerland),
 *     INRIA Sophia-Antipolis (France),
 *     Max-Planck-Institute Saarbruecken (Germany),
 *     and Tel-Aviv University (Israel).  All rights reserved.
 *
 *     Permission is hereby granted, free of charge, to any person obtaining a
 *     copy of this software and associated documentation files (the
 *     "Software"), to deal in the Software without restriction, including
 *     without limitation the rights to use, copy, modify, merge, publish,
 *     distribute, sublicense, and/or sell copies of the Software, and to
 *     permit persons to whom the Software is furnished to do so, subject to
 *     the following conditions:
 *
 *     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 *     OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 *     CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 *     TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * */

#ifndef POLYHEDRON_DEMO_STL_PLUGIN_IMP_H
#define POLYHEDRON_DEMO_STL_PLUGIN_IMP_H

// Copyright see CGAL-4.3
// Amendments by PTW (NUMA)
#include <CGAL/Modifier_base.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <fstream>

#include "mesher_cgal.h"
#include <boost/tuple/tuple.hpp>

template <class Kernel, class Polyhedron>
void Build_from_stl<Kernel, Polyhedron>::operator()( HDS& hds) {
      if(!read(is, meshPoints, this->mesh)) return;

      CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
      B.begin_surface( meshPoints.size(), this->mesh.size());
      typedef typename Points_3::size_type size_type;
      
      for(size_type i=0; i < meshPoints.size(); i++){
        B.add_vertex( meshPoints[i]);
      }
      for(size_type i=0; i < this->mesh.size(); i++){
        B.begin_facet(); 
        B.add_vertex_to_facet( this->mesh[i].template get<0>());
        B.add_vertex_to_facet( this->mesh[i].template get<1>());
        B.add_vertex_to_facet( this->mesh[i].template get<2>());
        B.end_facet();
      }
      if(B.error())
        {
          std::cerr << "An error occured while creating a Polyhedron" << std::endl;
          B.rollback();
        }
      
      B.end_surface();
}

template <class Kernel, class Polyhedron>
bool Build_from_stl<Kernel, Polyhedron>::read(std::istream& input, Points_3& points, Surface& surface, int /*offset*/) {
      std::string s, solid("solid"), facet("facet"), outer("outer"), loop("loop"), vertex("vertex"), endloop("endloop"), endsolid("endsolid");

      std::map<Point_3, int> vertex_index;
      int index = 0;
      int ijk[3];
      Point_3 p;

      input >> s;
      if(s == solid){
        std::getline(input, s);
      } else {
        std::cerr << "We expect keyword 'solid'" << std::endl;
        return false;
      }

      while(input >> s){
        if(s == endsolid){
          //std::cerr << "found endsolid" << std::endl;
        } else if(s == facet){
          //std::cerr << "found facet" << std::endl;
          std::getline(input, s); // ignore the normal
          input >> s;
          if(s != outer){
            std::cerr << "Expect 'outer' and got " << s << std::endl;
            return false;
          }
          input >> s;
          if(s != loop){
            std::cerr << "Expect 'loop' and got " << s << std::endl;
            return false;
         }
          int count = 0;
          do {
            input >> s;
            if(s == vertex){
              //      std::cerr << "found vertex" << std::endl;
              if(count < 3){
                input >> p;
                if(vertex_index.find(p) == vertex_index.end()){
                  ijk[count] = index;
                  vertex_index[p] = index++;
                  points.push_back(p);
                } else {
                  ijk[count] = vertex_index[p];
                }
                ++count;
              } else {
                std::cerr << "We can only read triangulated surfaces" << std::endl;
                return false;
              }
            }
          }while(s != endloop);

          surface.push_back(boost::make_tuple(ijk[0], ijk[1], ijk[2]));
        }
      }
      return true;
}

template <class Kernel, class Polyhedron>
int Build_from_stl<Kernel, Polyhedron>::load_stl(Polyhedron& P) {

      // Try to read STL file in a polyhedron
      P.delegate(*this);
     
      if(!P.is_valid()){
        std::cerr << "Error: Invalid polyhedron" << std::endl;
        return -2;
      }  

      if(P.empty()){
        std::cerr << "Error: Polyhedron empty" << std::endl;
        return -1;
      }  

      return 0;
}

#endif
