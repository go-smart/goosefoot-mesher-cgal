// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stéphane Tayeb (original Labeled_mesh_domain_3.h)
// Author(s)     : Phil Weir (NUMA Engineering modifications)
//
//******************************************************************************
// File Description :
// class Signed_mesh_domain_3. See class description.
//******************************************************************************


#ifndef CGAL_MESH_3_SIGNED_MESH_DOMAIN_3_H
#define CGAL_MESH_3_SIGNED_MESH_DOMAIN_3_H

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4180) // qualifier applied to function type has no meaning; ignored
#endif

#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>

namespace CGAL {


/**
 * \class Signed_mesh_domain_3
 *
 * The range of function f is the integers, Z.
 * Let p be a Point.
 *  - f(p)<=0 means that p is outside domain.
 *  - f(p) =a, a!=0 means that p is inside subdomain a.
 *
 *  Any boundary facet is labelled <a,b>, a<b, where a and b are the
 *  tags of its incident subdomain.
 *  Thus, a boundary facet of the domain is labelled <a,b>, where a<=0, b!=0.
 */
template<class Function, class BGT>
class Signed_mesh_domain_3
 : public Mesh_3::Labeled_mesh_domain_3<Function, BGT >
{
public:
  /// Base type
  typedef Mesh_3::Labeled_mesh_domain_3<Function, BGT> Base;

  /// Public types
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Subdomain Subdomain;
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Subdomain_index Subdomain_index;
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Surface_patch Surface_patch;
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Segment_3  Segment_3;
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Ray_3      Ray_3;
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Line_3     Line_3;
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Vector_3     Vector_3;
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Surface_patch_index     Surface_patch_index;
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Index     Index;
  typedef typename Mesh_3::Labeled_mesh_domain_3<Function, BGT>::Intersection     Intersection;
  typedef typename BGT::Iso_cuboid_3  Iso_cuboid_3;
  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::Point_3 Point_3;
  typedef typename Base::FT FT;

  /**
   * @brief Constructor
   */
  Signed_mesh_domain_3(const Function& f,
                         const Sphere_3& bounding_sphere,
                         const FT& error_bound = FT(1e-3),
                         CGAL::Random* p_rng = NULL);

  Signed_mesh_domain_3(const Function& f,
                         const Bbox_3& bbox,
                         const FT& error_bound = FT(1e-3),
                         CGAL::Random* p_rng = NULL);

  Signed_mesh_domain_3(const Function& f,
                         const Iso_cuboid_3& bbox,
                         const FT& error_bound = FT(1e-3),
                         CGAL::Random* p_rng = NULL);

  /// Destructor
  virtual ~Signed_mesh_domain_3()
  {
    if(delete_rng_)
      delete p_rng_;
  }

  /**
   * Constructs  a set of \ccc{n} points on the surface, and output them to
   *  the output iterator \ccc{pts} whose value type is required to be
   *  \ccc{std::pair<Points_3, Index>}.
   */
  struct Construct_initial_points
  {
    Construct_initial_points(const Signed_mesh_domain_3& domain)
      : r_domain_(domain) {}

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 12) const;

  private:
    const Signed_mesh_domain_3& r_domain_;
  };

  /// Returns Construct_initial_points object
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  struct Is_in_domain
  {
    Is_in_domain(const Signed_mesh_domain_3& domain) : r_domain_(domain) {}

    Subdomain operator()(const Point_3& p) const
    {
      // f(p)<=0 means p is outside the domain
      Subdomain_index index = (r_domain_.signed_function_)(p);
      if ( Subdomain_index() >= index )
        return Subdomain();
      else
        return Subdomain(index);
    }
  private:
    const Signed_mesh_domain_3& r_domain_;
  };

  Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }

  /**
   * Returns true is the element \ccc{type} intersect properly any of the
   * surface patches describing the either the domain boundary or some
   * subdomain boundary.
   * \ccc{Type} is either \ccc{Segment_3}, \ccc{Ray_3} or \ccc{Line_3}.
   * Parameter index is set to the index of the intersected surface patch
   * if \ccc{true} is returned and to the default \ccc{Surface_patch_index}
   * value otherwise.
   */
  struct Do_intersect_surface
  {
    Do_intersect_surface(const Signed_mesh_domain_3& domain)
      : r_domain_(domain) {}

    Surface_patch operator()(const Segment_3& s) const
    {
      return this->operator()(s.source(), s.target());
    }

    Surface_patch operator()(const Ray_3& r) const
    {
      return clip_to_segment(r);
    }

    Surface_patch operator()(const Line_3& l) const
    {
      return clip_to_segment(l);
    }

  private:
    /// Returns true if points \c a & \c b do not belong to the same subdomain
    /// \c index is set to the surface index of subdomains f(a), f(b)
    Surface_patch operator()(const Point_3& a, const Point_3& b) const
    {
      // If f(a) != f(b), then [a,b] intersects some surface. Here we consider
      // [a,b] intersects surface_patch labelled <f(a),f(b)> (or <f(b),f(a)>).
      // It may be false, further rafinement will improve precision
      const Subdomain_index value_a = r_domain_.signed_function_(a);
      const Subdomain_index value_b = r_domain_.signed_function_(b);

      if ( (value_a > 0 || value_b > 0) && value_a != value_b )
        return Surface_patch(r_domain_.make_surface_index(value_a, value_b));
      else
        return Surface_patch();
    }

    /**
     * Clips \c query to a segment \c s, and call operator()(s)
     */
    template<typename Query>
    Surface_patch clip_to_segment(const Query& query) const
    {
      typename cpp11::result_of<typename BGT::Intersect_3(Query, Iso_cuboid_3)>::type
        clipped = CGAL::intersection(query, r_domain_.bounding_box());

      if(clipped)
#if CGAL_INTERSECTION_VERSION > 1
        if(const Segment_3* s = boost::get<Segment_3>(&*clipped))
          return this->operator()(*s);
#else
        if(const Segment_3* s = object_cast<Segment_3>(&clipped))
          return this->operator()(*s);
#endif
        
      return Surface_patch();
    }

  private:
    const Signed_mesh_domain_3& r_domain_;
  };

  /// Returns Do_intersect_surface object
  Do_intersect_surface do_intersect_surface_object() const
  {
    return Do_intersect_surface(*this);
  }

  /**
   * Returns a point in the intersection of the primitive \ccc{type}
   * with some boundary surface.
   * \ccc{Type1} is either \ccc{Segment_3}, \ccc{Ray_3} or \ccc{Line_3}.
   * The integer \ccc{dimension} is set to the dimension of the lowest
   * dimensional face in the input complex containing the returned point, and
   * \ccc{index} is set to the index to be stored at a mesh vertex lying
   * on this face.
   */
  struct Construct_intersection
  {
    Construct_intersection(const Signed_mesh_domain_3& domain)
      : r_domain_(domain) {}

    Intersection operator()(const Segment_3& s) const
    {
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      CGAL_precondition(r_domain_.do_intersect_surface_object()(s));
#endif // NOT CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      return this->operator()(s.source(),s.target());
    }

    Intersection operator()(const Ray_3& r) const
    {
      return clip_to_segment(r);
    }

    Intersection operator()(const Line_3& l) const
    {
      return clip_to_segment(l);
    }

  private:
    /**
     * Returns a point in the intersection of [a,b] with the surface
     * \c a must be the source point, and \c b the out point. It's important
     * because it drives bisection cuts.
     * Indeed, the returned point is the first intersection from \c [a,b]
     * with a subdomain surface.
     */
    Intersection operator()(const Point_3& a, const Point_3& b) const
    {
      // Functors
      typename BGT::Compute_squared_distance_3 squared_distance =
                                      BGT().compute_squared_distance_3_object();
      typename BGT::Construct_midpoint_3 midpoint =
                                      BGT().construct_midpoint_3_object();

      // Non const points
      Point_3 p1 = a;
      Point_3 p2 = b;
      Point_3 mid = midpoint(p1, p2);

      // Cannot be const: those values are modified below.
      Subdomain_index value_at_p1 = r_domain_.signed_function_(p1);
      Subdomain_index value_at_p2 = r_domain_.signed_function_(p2);
      Subdomain_index value_at_mid = r_domain_.signed_function_(mid,true);

      // If both extremities are in the same subdomain,
      // there is no intersection.
      // This should not happen...
      if( value_at_p1 == value_at_p2 || (value_at_p1 <= 0 && value_at_p2 <= 0) )
      {
        return Intersection();
      }

      // Construct the surface patch index and index from the values at 'a'
      // and 'b'. Even if the bissection find out a different pair of
      // values, the reported index will be constructed from the initial
      // values.
      const Surface_patch_index sp_index =
        r_domain_.make_surface_index(value_at_p1, value_at_p2);
      const Index index = r_domain_.index_from_surface_patch_index(sp_index);

      // Else lets find a point (by bisection)
      // Bisection ends when the point is near than error bound from surface
      while(true)
      {
        // If the two points are enough close, then we return midpoint
        if ( squared_distance(p1, p2) < r_domain_.squared_error_bound_ )
        {
          CGAL_assertion(value_at_p1 != value_at_p2 && (value_at_p1 > 0 || value_at_p2 > 0));
          return Intersection(mid, index, 2);
        }

        // Else we must go on
        // Here we consider that p1(a) is the source point. Thus, we keep p1 and
        // change p2 if f(p1)!=f(p2).
        // That allows us to find the first intersection from a of [a,b] with
        // a surface.
        if ( value_at_p1 != value_at_mid && (value_at_p1 > 0 || value_at_mid > 0) )
        {
          p2 = mid;
          value_at_p2 = value_at_mid;
        }
        else
        {
          p1 = mid;
          value_at_p1 = value_at_mid;
        }

        mid = midpoint(p1, p2);
        value_at_mid = r_domain_.signed_function_(mid,true);
      }
    }

    /// Clips \c query to a segment \c s, and call operator()(s)
    template<typename Query>
    Intersection clip_to_segment(const Query& query) const
    {
      typename cpp11::result_of<typename BGT::Intersect_3(Query, Iso_cuboid_3)>::type
        clipped = CGAL::intersection(query, r_domain_.bounding_box());

      if(clipped)
#if CGAL_INTERSECTION_VERSION > 1
        if(const Segment_3* s = boost::get<Segment_3>(&*clipped))
          return this->operator()(*s);
#else
        if(const Segment_3* s = object_cast<Segment_3>(&clipped))
          return this->operator()(*s);
#endif
      
      return Intersection();
    }

  private:
    const Signed_mesh_domain_3& r_domain_;
  };

  /// Returns Construct_intersection object
  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }

  /// Returns the bounding sphere of an Iso_cuboid_3
  Sphere_3 bounding_sphere(const Iso_cuboid_3& bbox) const
  {
    typename BGT::Construct_sphere_3 sphere = BGT().construct_sphere_3_object();
    return sphere((bbox.min)(), (bbox.max)());
  }

protected:
  /// The function which answers subdomain queries (this would not be required if Lmd3._function were protected not private)
  const Function signed_function_;
  FT squared_error_bound_;

private:
  CGAL::Random* p_rng_;
  bool delete_rng_;

private:
  // Disabled copy constructor & assignment operator
  typedef Signed_mesh_domain_3<Function,BGT> Self;
  Signed_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

private:
  /// Returns Surface_patch_index from \c i and \c j
  Surface_patch_index make_surface_index(const Subdomain_index i,
                                   const Subdomain_index j) const
  {
    if ( i < j ) return Surface_patch_index(i,j);
    else return Surface_patch_index(j,i);
  }

  /// Returns squared error bound from \c sphere and \c error
  FT squared_error_bound(const Sphere_3& sphere, const FT& error) const
  {
    typename BGT::Compute_squared_radius_3 squared_radius =
                                    BGT().compute_squared_radius_3_object();
    return squared_radius(sphere)*error*error;
  }

};  // end class Signed_mesh_domain_3


//-------------------------------------------------------
// Method implementation
//-------------------------------------------------------
template<class F, class BGT>
Signed_mesh_domain_3<F,BGT>::Signed_mesh_domain_3(
                       const F& f,
                       const Sphere_3& bounding_sphere,
                       const FT& error_bound,
                       CGAL::Random* p_rng)
: CGAL::Mesh_3::Labeled_mesh_domain_3<F, BGT>(f, bounding_sphere, error_bound)
, signed_function_(f)
, squared_error_bound_(squared_error_bound(bounding_sphere,error_bound))
, p_rng_(p_rng)
, delete_rng_(false)
{
  // TODO : CGAL_ASSERT(0 < f(bounding_sphere.get_center()) ) ?
  if(!p_rng_)
  {
    p_rng_ = new CGAL::Random(0);
    delete_rng_ = true;
  }
}

template<class F, class BGT>
Signed_mesh_domain_3<F,BGT>::Signed_mesh_domain_3(
                       const F& f,
                       const Bbox_3& bbox,
                       const FT& error_bound,
                       CGAL::Random* p_rng)
: CGAL::Mesh_3::Labeled_mesh_domain_3<F, BGT>(f, bbox, error_bound)
, signed_function_(f)
, squared_error_bound_(squared_error_bound(CGAL::Mesh_3::Labeled_mesh_domain_3<F, BGT>::bounding_box(),error_bound))
, p_rng_(p_rng)
, delete_rng_(false)
{
  // TODO : CGAL_ASSERT(0 < f(bounding_sphere.get_center()) ) ?
  if(!p_rng_)
  {
    p_rng_ = new CGAL::Random(0);
    delete_rng_ = true;
  }
}

template<class F, class BGT>
Signed_mesh_domain_3<F,BGT>::Signed_mesh_domain_3(
                       const F& f,
                       const Iso_cuboid_3& bbox,
                       const FT& error_bound,
                       CGAL::Random* p_rng)
: CGAL::Mesh_3::Labeled_mesh_domain_3<F, BGT>(f, bbox, error_bound, p_rng)
, signed_function_(f)
, squared_error_bound_(squared_error_bound(CGAL::Mesh_3::Labeled_mesh_domain_3<F, BGT>::bounding_box(),error_bound))
, p_rng_(p_rng)
, delete_rng_(false)
{
  // TODO : CGAL_ASSERT(0 < f( bbox.get_center()) ) ?
  if(!p_rng_)
  {
    p_rng_ = new CGAL::Random(0);
    delete_rng_ = true;
  }
}

template<class F, class BGT>
template<class OutputIterator>
OutputIterator
Signed_mesh_domain_3<F,BGT>::Construct_initial_points::operator()(
                                                    OutputIterator pts,
                                                    const int nb_points) const
{
  // Create point_iterator on and in bounding_sphere
  typedef Random_points_on_sphere_3<Point_3> Random_points_on_sphere_3;
  typedef Random_points_in_sphere_3<Point_3> Random_points_in_sphere_3;


  const FT squared_radius = BGT().compute_squared_radius_3_object()(
      r_domain_.bounding_sphere(r_domain_.bounding_box()));

  const double radius = std::sqrt(CGAL::to_double(squared_radius));

  CGAL::Random& rng = *(r_domain_.p_rng_);
  Random_points_on_sphere_3 random_point_on_sphere(radius, rng);
  Random_points_in_sphere_3 random_point_in_sphere(radius, rng);

  // Get some functors
  typename BGT::Construct_segment_3 segment_3 =
                              BGT().construct_segment_3_object();
  typename BGT::Construct_vector_3 vector_3 =
                              BGT().construct_vector_3_object();
  typename BGT::Construct_translated_point_3 translate =
                              BGT().construct_translated_point_3_object();
  typename BGT::Construct_center_3 center = BGT().construct_center_3_object();

  // Get translation from origin to sphere center
  Point_3 center_pt = center(r_domain_.bounding_sphere(r_domain_.bounding_box()));
  const Vector_3 sphere_translation = vector_3(CGAL::ORIGIN, center_pt);

  // Create nb_point points
  int n = nb_points;
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "construct initial points:\n";
#endif
  while ( 0 != n )
  {
    // Get a random segment
    const Point_3 random_point = translate(*random_point_on_sphere,
                                           sphere_translation);
    const Segment_3 random_seg = segment_3(center_pt, random_point);

    // Add the intersection to the output if it exists
    Surface_patch surface = r_domain_.do_intersect_surface_object()(random_seg);
    if ( surface )
    {
      const Point_3 intersect_pt = CGAL::cpp11::get<0>(
          r_domain_.construct_intersection_object()(random_seg));
      *pts++ = std::make_pair(intersect_pt,
                              r_domain_.index_from_surface_patch_index(*surface));
      --n;

#ifdef CGAL_MESH_3_VERBOSE
      std::cerr << boost::format("\r             \r"
                                 "%1%/%2% initial point(s) found...")
                   % (nb_points - n)
                   % nb_points;
#endif
    }
    else
    {
      // Get a new random point into sphere as center of object
      // It may be necessary if the center of the domain is empty, e.g. torus
      // In general case, it is good for input point dispersion
      ++random_point_in_sphere;
      center_pt = translate(*random_point_in_sphere, sphere_translation);
    }
    ++random_point_on_sphere;
  }

#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "\n";
#endif
  return pts;
}


}  // end namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_SIGNED_MESH_DOMAIN_3_H
