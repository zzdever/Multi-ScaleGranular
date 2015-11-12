#ifndef Ray_h
#define Ray_h

#include "Vector3.h"

namespace FastBVH{
	
struct Ray {
  FastBVH::Vector3 o; // Ray Origin
  FastBVH::Vector3 d; // Ray Direction
  FastBVH::Vector3 inv_d; // Inverse of each Ray Direction component

  Ray(const FastBVH::Vector3& o, const FastBVH::Vector3& d)
    : o(o), d(d), inv_d(FastBVH::Vector3(1,1,1).cdiv(d)) { }
};

} //namespace BVH

#endif
