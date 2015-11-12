#ifndef IntersectionInfo_h_
#define IntersectionInfo_h_

#include "Object.h"

namespace FastBVH{
	
struct Object;

struct IntersectionInfo {
  float t; // Intersection distance along the ray
  const FastBVH::Object* object; // Object that was hit
  FastBVH::Vector3 hit; // Location of the intersection
};

} //namespace FastBVH

#endif
