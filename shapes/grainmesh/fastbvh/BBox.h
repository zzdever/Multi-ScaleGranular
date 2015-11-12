#ifndef BBox_h
#define BBox_h

#include "Ray.h"
#include "Vector3.h"
#include <stdint.h>

struct BBox {
  FastBVH::Vector3 min, max, extent;
  BBox() { }
  BBox(const FastBVH::Vector3& min, const FastBVH::Vector3& max);
  BBox(const FastBVH::Vector3& p);

  bool intersect(const FastBVH::Ray& ray, float *tnear, float *tfar) const;
  void expandToInclude(const FastBVH::Vector3& p);
  void expandToInclude(const BBox& b);
  uint32_t maxDimension() const;
  float surfaceArea() const;
};

#endif
