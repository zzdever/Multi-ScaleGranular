
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
                m_spherePack[i].d/2)
			);
            // objects may be modified
            objects_bak.push_back(objects.at(i));
		}
		
		// Compute a BVH for this object set
        bvh = new BVH(&objects);
	}

    SpherePack::Sphere* rayIntersect(Point o, Point d, FastBVH::IntersectionInfo &I){
        FastBVH::Ray ray(FastBVH::Vector3(o.x, o.y, o.z), normalize(FastBVH::Vector3(d.x, d.y, d.z)));
        if (! bvh->getIntersection(ray, &I, false))
            return NULL;

        if(I.t < 0)
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
