#ifndef BVH_h
#define BVH_h

#include "BBox.h"
#include <vector>
#include <stdint.h>
#include "Object.h"
#include "IntersectionInfo.h"
#include "Ray.h"

//! Node descriptor for the flattened tree
struct BVHFlatNode {
  BBox bbox;
  uint32_t start, nPrims, rightOffset;
};

//! \author Brandon Pelfrey
//! A Bounding Volume Hierarchy system for fast Ray-Object intersection tests
class BVH {
  uint32_t nNodes, nLeafs, leafSize;
  std::vector<FastBVH::Object*>* build_prims;

  //! Build the BVH tree out of build_prims
  void build();

  // Fast Traversal System
  BVHFlatNode *flatTree;

  public:
  BVH(std::vector<FastBVH::Object*>* objects, uint32_t leafSize=4);
  bool getIntersection(const FastBVH::Ray& ray, FastBVH::IntersectionInfo *intersection, bool occlusion) const ;

  ~BVH();
};

#endif
