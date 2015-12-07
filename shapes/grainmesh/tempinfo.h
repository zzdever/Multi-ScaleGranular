
#ifndef TEMPINFO_H
#define TEMPINFO_H

#include <mitsuba/core/plugin.h>
#include "bvhinterface.h"


MTS_NAMESPACE_BEGIN

typedef enum {
    EPT,
    VPT,
    DA,
    SCATTEROMETER,
    DEBUG
} MODES_T;

struct SphereIntersectionInfo{
    FastBVH::IntersectionInfo I;
    SpherePack::Sphere* sphere;
    int index;
};

struct GrainIntersectionInfo{
    MODES_T rendermode;

    const Shape* shape;

    Point2 uv;
    Normal n;
    Vector dpdu, dpdv;
};

struct ScatterometerReturnInfo{
    bool isIntersected;
    Point intersection;
    Vector outDirection;
};



struct hitInfo{
    void * hitObject;
    Point hitPoint;
    Float u,v;
};


MTS_NAMESPACE_END

#endif
