#include <mitsuba/core/plugin.h>


MTS_NAMESPACE_BEGIN

struct IntersectedSphereInfo{
    Float radius;
    Float x;
    Float y;
    Float z;
    Float w;
    Point center;
    int name;
    Float ox, oy, oz, t;
    Float ix, iy, iz;
public:
    IntersectedSphereInfo():radius(0), x(0), y(0), z(0), w(0) {}
};

MTS_NAMESPACE_END
