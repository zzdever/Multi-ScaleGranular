#ifndef Sphere_h_
#define Sphere_h_

#include <cmath>
#include "Object.h"

namespace FastBVH{

struct Sphere : public FastBVH::Object {
  FastBVH::Vector3 center; // Center of the sphere
  float r, r2; // Radius, Radius^2

  Sphere(const FastBVH::Vector3& center, float radius)
    : center(center), r(radius), r2(radius*radius) { }

  bool getIntersection(const FastBVH::Ray& ray, FastBVH::IntersectionInfo* I) const {
    FastBVH::Vector3 s = center - ray.o;
    float sd = s * ray.d;
    float ss = s * s;

    // Compute discriminant
    float disc = sd*sd - ss + r2;

    // Complex values: No intersection
    if( disc < 0.f ) return false;

    // Assume we are not in a sphere... The first hit is the lesser valued
    I->object = this;
    I->t = sd - sqrt(disc);
    return true;
  }

  FastBVH::Vector3 getNormal(const FastBVH::IntersectionInfo& I) const {
    return normalize(I.hit - center);
  }

  BBox getBBox() const {
    return BBox(center-FastBVH::Vector3(r,r,r), center+FastBVH::Vector3(r,r,r));
  }

  FastBVH::Vector3 getCentroid() const {
    return center;
  }

};

} // namespace FastBVH

#endif
