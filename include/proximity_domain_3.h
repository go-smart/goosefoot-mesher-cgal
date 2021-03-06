#ifndef PROXIMITY_DOMAIN_3_H
#define PROXIMITY_DOMAIN_3_H

#include "mesher_cgal.h"

#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef K::Iso_cuboid_3 Iso_cuboid;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

namespace mesherCGAL {
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
      Proximity_domain_field_3(Tree* tree, double near_cl, double far_cl, std::vector< std::unique_ptr<Zone> >& zones, float zone_radius, float centre_radius,
                    const Point* centre, Tree* needle_tree, float granularity, const Iso_cuboid& bbox, bool solid_zones) :
                            max_(CGAL::Max<FT>()), min_(CGAL::Min<FT>()),
                            tree_(tree), near_cl_(near_cl), far_cl_(far_cl), zones_(zones), zone_radius_(zone_radius),
                            centre_radius_(centre_radius),
                            centre_(centre), needle_tree_(needle_tree), granularity_(granularity),
                            bbox_(bbox), checks_(0), total_checks_(0), values_size_(0), solid_zones_(solid_zones)
        {
            max_dist[0] = 0.; max_dist[1] = 0.; max_dist[2] = 0.;

            minx = (int)(bbox.xmin()/granularity_); nx = 1 + (int)(bbox.xmax()/granularity_) - minx;
            miny = (int)(bbox.ymin()/granularity_); ny = 1 + (int)(bbox.ymax()/granularity_) - miny;
            minz = (int)(bbox.zmin()/granularity_); nz = 1 + (int)(bbox.zmax()/granularity_) - minz;
            length_scale = std::pow(bbox.volume(), 1./3);

            std::cout << "Bounding Box coordinates: "
                << "[" << bbox.xmin() << ", " << bbox.ymin() << ", " << bbox.zmin() << "] -> "
                << "[" << bbox.xmax() << ", " << bbox.ymax() << ", " << bbox.zmax() << "]"
                << std::endl;
            std::cout << "Length scale: " << length_scale << std::endl;
            std::cout << "CL bitmap: "
                << "[" << minx << ", " << miny << ", " << minz << "] + "
                << "[" << nx << ", " << ny << ", " << nz << "]"
                << std::endl;
            std::cout << "Array size: " << nx * ny * nz << std::endl;

            //visualization_create_structured_grid(minx, miny, minz, nx, ny, nz, granularity_);

            values_ = (FT*)calloc(sizeof(FT), nx*ny*nz);
            dists_ = (FT*)calloc(sizeof(FT), nx*ny*nz);
            max_cl_ = -1.;
        }
      //RMV : insert destructor to clear values_ allocations

      /* http://www.vtk.org/Wiki/VTK/Examples/Cxx/ImageData/IterateImageData */
      void output_characteristic_lengths(std::string filename)
      {
          vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
          imageData->SetDimensions(nx, ny, nz);
          imageData->SetNumberOfScalarComponents(2);
          imageData->SetScalarTypeToDouble();
          imageData->SetSpacing(granularity_, granularity_, granularity_);
          imageData->SetOrigin(minx * granularity_, miny * granularity_, minz * granularity_);

          for (int z = 0 ; z < nz ; z++)
              for (int y = 0 ; y < ny ; y++)
                  for (int x = 0 ; x < nx ; x++) {
                      double *pix = static_cast<double*>(imageData->GetScalarPointer(x, y, z));
                      int ix = nx*(ny*z + y) + x;
                      pix[0] = values_[ix] > 1e-20 ? values_[ix] : -1.0;
                      pix[1] = dists_[ix] > 1e-20 ? dists_[ix] : -2.0;
                  }

          vtkSmartPointer<vtkXMLImageDataWriter> imageWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
          imageWriter->SetInput(imageData);
          imageWriter->SetFileName(filename.c_str());
          imageWriter->Update();
          imageWriter->Write();
      }

      FT to_cl(FT closest_cl, FT dist) const
      {
          FT scaling = 1;

          if (dist <= 4 * closest_cl * closest_cl + 1e-10)
              scaling = 0;
          else
              scaling *= pow((dist - closest_cl * closest_cl * 4)/(length_scale * length_scale * 0.25), 0.5);

          return (far_cl_ - (far_cl_-closest_cl)*(exp(-scaling) - exp(-1)) / (1 - exp(-1)));
      }

      /// Returns size
      FT operator()(const Point_3& q, const int , const Index& ) const
      {
          Point_3 p(
                  min_(max_(q.x(), bbox_.xmin()), bbox_.xmax()),
                  min_(max_(q.y(), bbox_.ymin()), bbox_.ymax()),
                  min_(max_(q.z(), bbox_.zmin()), bbox_.zmax())
          );
          int x = (int)(p.x()/granularity_);
          int y = (int)(p.y()/granularity_);
          int z = (int)(p.z()/granularity_);
          //int ix = nx*(ny*(min_(max_(z - minz, 0), nz - 1)) + min_(max_(y - miny, 0), ny - 1)) + min_(max_(x - minx, 0), nx - 1);
          int ix = nx*(ny*(z - minz) + y - miny) + x - minx;

          FT cl = values_[ix];

          checks_++;

          if ( checks_ % 10000 == 0 ) {
              total_checks_ += checks_;
              std::cout << "Cached " << values_size_
                  << " blocks (miss rate " << values_size_*100./total_checks_
                  << "%) - max. cl: " << max_cl_
                  << " - max. distances from key points/surfaces: " << max_dist[0]
                  << " " << max_dist[1] << " " << max_dist[2] << std::endl;
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

          /* Everywhere dist appears, it is squared */
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

          for (auto&& zone : zones_) {
              FT zone_dist = -1;
              if (solid_zones_ && zone->is_container() && zone->contains(p)) {
                  zone_dist = 0;
              }
              else if (zone->has_tree()) {
                  zone_dist = zone->squared_distance(p);
                  if (zone_dist > zone_radius_ * zone_radius_)
                      zone_dist = zone_dist + zone_radius_ * zone_radius_ - 2. * sqrt(zone_dist) * zone_radius_;
                  else zone_dist = 0.;
                  if (zone_dist < 0) zone_dist = 0.;
              }
              //zone_dist = CGAL::min(zone_dist - near_cl_ * 2, 0.0); /* Introduce near_cl region around zone */
              if (zone_dist >= 0 && (dist < 0 || (to_cl(closest_cl, dist) > to_cl(zone->get_cl(), zone_dist)))) {
                  dist = zone_dist;
                  closest_cl = zone->get_cl();
              }
          }

          if (centre_ != NULL) {
              FT centre_dist = squared_distance(*centre_, p);
              if (centre_dist > centre_radius_ * centre_radius_)
                  centre_dist = centre_dist + centre_radius_ * centre_radius_ - 2. * sqrt(centre_dist) * centre_radius_;
              else centre_dist = 0.;
              if (centre_dist < 0) centre_dist = 0.;
              if (dist < 0 || to_cl(closest_cl, dist) > to_cl(near_cl_, centre_dist)) {
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
          dists_[ix] = dist;
          //visualization_set_allocation_order(ix, values_size_, cl, needle_dist, p.x(), p.y(), p.z());

          return cl;
      }

      mutable float max_cl_;//RMV
    private:
      Tree* tree_;
      mutable FT* values_;
      mutable FT* dists_;
      FT near_cl_, far_cl_;
      std::vector< std::unique_ptr<Zone> >& zones_;
      FT zone_radius_, centre_radius_;
      const Point* centre_;
      Tree* needle_tree_;
      FT granularity_;
      const Iso_cuboid& bbox_;
      mutable int checks_, total_checks_, values_size_;
      mutable FT max_dist[3];
      bool solid_zones_;
    };
}

#endif
