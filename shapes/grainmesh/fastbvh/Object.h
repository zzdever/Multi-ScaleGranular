#ifndef Object_h_
#define Object_h_

#include "IntersectionInfo.h"
#include "Ray.h"
#include "BBox.h"

namespace FastBVH{

struct Object {
  //! All "Objects" must be able to test for intersections with rays.
  virtual bool getIntersection(
      const FastBVH::Ray& ray,
      FastBVH::IntersectionInfo* intersection)
    const = 0;

  //! Return an object normal based on an intersection
  virtual FastBVH::Vector3 getNormal(const FastBVH::IntersectionInfo& I) const = 0;

  //! Return a bounding box for this object
  virtual BBox getBBox() const = 0;

  //! Return the centroid for this object. (Used in BVH Sorting)
  virtual FastBVH::Vector3 getCentroid() const = 0;
};

} //namespace BVH

#endif
