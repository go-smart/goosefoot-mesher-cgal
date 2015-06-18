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
  typedef typename Base::Sphere_3 Sphere_3;
  typedef typename Base::Point_3 Point_3;
  typedef typename Base::FT FT;

  /**
   * Constructor
   * @param f the function which negative values defines the domain
   * @param bounding_sphere a bounding sphere of the domain
   * @param error_bound the error bound relative to the sphere radius
   */
  Signed_mesh_domain_3(const Function& f,
                         const Sphere_3& bounding_sphere,
                         const FT& error_bound = FT(1e-3))
    : Base(f, bounding_sphere, error_bound), signed_function_(f)  {}

  /// Destructor
  virtual ~Signed_mesh_domain_3() {}

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

  //TODO: Does the interval bisection get affected by artificial subdomains outside the domain?

protected:
  /// The function which answers subdomain queries (this would not be required if Lmd3._function were protected not private)
  const Function signed_function_;

private:
  // Disabled copy constructor & assignment operator
  typedef Signed_mesh_domain_3<Function,BGT> Self;
  Signed_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Signed_mesh_domain_3


}  // end namespace CGAL

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif // CGAL_SIGNED_MESH_DOMAIN_3_H
