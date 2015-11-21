
#include <vector>

#include "fastbvh/BVH.h"
#include "fastbvh/Sphere.h"
#include "fastbvh/IntersectionInfo.h"

#include "sphere_pack.h"

using std::vector;
using namespace FastBVH;

class BVHInterface{
public:
    BVHInterface(struct SpherePack::Sphere* pack, int count) : m_spherePack(pack), N(count){
        this->build();
    }

	inline void build(){
        // N = sizeof(m_spherePack) / sizeof(struct SpherePack::Sphere);
        for(size_t i=0; i<N; ++i) {
            objects.push_back(new FastBVH::Sphere(
                FastBVH::Vector3(m_spherePack[i].x, m_spherePack[i].y, m_spherePack[i].z),
                m_spherePack[i].d/2.0)
			);
            // objects may be modified
            objects_bak.push_back(objects.at(i));
		}
		
		// Compute a BVH for this object set
        bvh = new BVH(&objects);
	}

    SpherePack::Sphere* rayIntersect(Point _o, Point _d, FastBVH::IntersectionInfo &I){
        FastBVH::Vector3 o = FastBVH::Vector3(_o.x, _o.y, _o.z);
        FastBVH::Vector3 d = normalize(FastBVH::Vector3(_d.x, _d.y, _d.z));
        FastBVH::Ray ray(o, d);
        if (! bvh->getIntersection(ray, &I, false))
            return NULL;

        // is the origin in a sphere?
        if(length(o - static_cast<const Sphere*>(I.object)->center) <= m_spherePack[0].d/2.0) {
            o = o + I.t * d;
        }

        // trace again
        ray.o = o;
        if (! bvh->getIntersection(ray, &I, false))
            return NULL;

        if(I.t <= 0)
            return NULL;

        //printf("fbvh:ooo: %f, %f, %f\n fbvh:ddd: %f, %f, %f\n fbvh:hit: %f, %f, %f\n",
               //o.x, o.y, o.z, d.x, d.y, d.z, I.hit.x, I.hit.y, I.hit.z);

        for(size_t i=0; i<N; ++i) {
            if(objects_bak.at(i) == I.object)
                return &m_spherePack[i];
        }

        return NULL;
	}
	
private:
    BVH *bvh;
    struct SpherePack::Sphere* m_spherePack;
    unsigned int N;
    vector<FastBVH::Object*> objects, objects_bak;
};


# if 0
int main(int argc, char **argv) {
	for(size_t i=0; i<width; ++i) {
		for(size_t j=0; j<height; ++j) {
			size_t index = 3*(width * j + i);

			// This is only valid for square aspect ratio images
			Ray ray(camera_position, normalize(u*camera_u + v*camera_v + fov*camera_dir));

			IntersectionInfo I;
			bool hit = bvh.getIntersection(ray, &I, false);
            Vector3 normal = I.object->getNormal(I);
        }
    }
}
#endif
