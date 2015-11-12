#include <cstdio>
#include <vector>
#include <cstdlib>
#include <iostream>
#include "BVH.h"
#include "Sphere.h"

using std::vector;
using namespace FastBVH;

#include "../sphere_pack.h"
#include "../sphere_pack_data.h"

int main(int argc, char **argv) {
	int m_dictCount = sizeof(spherePacks) / sizeof(spherePacks[0]);
	int m_packCount = sizeof(spherePacks[0]) / sizeof(struct SpherePack::Sphere);
	std::cout<<m_dictCount<<" "<<m_packCount<<std::endl;
			
	struct SpherePack::Sphere* spherePack = (struct SpherePack::Sphere*)spherePacks[0];
	
	const unsigned int N = 80;
	vector<FastBVH::Object*> objects, objects_bak;
	printf("Constructing %d spheres...\n", N);
	for(size_t i=0; i<N; ++i) {
		objects.push_back(new FastBVH::Sphere(FastBVH::Vector3(spherePack[i].x, spherePack[i].y, spherePack[i].z), spherePack[i].d/2));
		objects_bak.push_back(objects.at(i));
	}

	// Compute a BVH for this object set
	BVH bvh(&objects);
	
	
	Vector3 center = Vector3(0.6762580508707393, 0.4620289160344975, 0.1241752283957101);
	
	Vector3 dir = Vector3(10,50,0.2285823979034143+0.000001);
	Vector3 origin = center - dir*100;
	//printf("%f,%f,%f\n", origin.x, origin.y, origin.z);
	//Ray ray(origin, normalize(Vector3(0,0,0) + dir));
	Ray ray(Vector3(-7.959488, 0.205154, 0.000048), normalize(Vector3(0.102441, 0.077679, 0.991701)));

	IntersectionInfo I;
	bool hit = bvh.getIntersection(ray, &I, false);
	printf("is hit: %d\n", hit);
	printf("hit point: %f,%f,%f\n", I.hit.x, I.hit.y, I.hit.z);
	const FastBVH::Sphere* s = static_cast<const FastBVH::Sphere*>(I.object);
	printf("canonical sphere center: %lf,%lf,%lf,%lf\n", s->center.x, s->center.y, s->center.z, s->r);
	
	
	struct SpherePack::Sphere* o;
	for(int i=0; i<N; ++i) {
		if(objects_bak.at(i) == I.object){
			std::cout<<i<<std::endl;
			printf("sphere center: %f,%f,%f\n", spherePack[i].x, spherePack[i].y, spherePack[i].z);			
			break;
		}
	   		//o = (struct SpherePack::Sphere*)&spherePacks[i];
	}
	
	//std::cout<<hit<<std::endl;
	
	
	
	return 0;
	
	

	// Allocate space for some image pixels
	const unsigned int width=800, height=800;
	float* pixels = new float[width*height*3];

	// Create a camera from position and focus point
	FastBVH::Vector3 camera_position(1.6, 1.3, 1.6);
	FastBVH::Vector3 camera_focus(0,0,0);
	FastBVH::Vector3 camera_up(0,1,0);

	// Camera tangent space
	FastBVH::Vector3 camera_dir = normalize(camera_focus - camera_position);
	FastBVH::Vector3 camera_u = normalize(camera_dir ^ camera_up);
	FastBVH::Vector3 camera_v = normalize(camera_u ^ camera_dir);

	printf("Rendering image (%dx%d)...\n", width, height);
	// Raytrace over every pixel
#pragma omp parallel for
	for(size_t i=0; i<width; ++i) {
		for(size_t j=0; j<height; ++j) {
			size_t index = 3*(width * j + i);

			float u = (i+.5f) / (float)(width-1) - .5f;
			float v = (height-1-j+.5f) / (float)(height-1) - .5f;
			float fov = .5f / tanf( 70.f * 3.14159265*.5f / 180.f);

			// This is only valid for square aspect ratio images
			FastBVH::Ray ray(camera_position, normalize(u*camera_u + v*camera_v + fov*camera_dir));

			IntersectionInfo I;
			bool hit = bvh.getIntersection(ray, &I, false);

			if(!hit) {
				pixels[index] = pixels[index+1] = pixels[index+2] = 0.f;
			} else {

				// Just for fun, we'll make the color based on the normal
				const FastBVH::Vector3 normal = I.object->getNormal(I);
				const FastBVH::Vector3 color(fabs(normal.x), fabs(normal.y), fabs(normal.z));

				pixels[index  ] = color.x;
				pixels[index+1] = color.y;
				pixels[index+2] = color.z;
			}
		}
	}

	// Output image file (PPM Format)
	printf("Writing out image file: \"render.ppm\"\n");
	FILE *image = fopen("render.ppm", "w");
	fprintf(image, "P6\n%d %d\n255\n", width, height);
	for(size_t j=0; j<height; ++j) {
		for(size_t i=0; i<width; ++i) {
			size_t index = 3*(width * j + i);
			unsigned char r = std::max(std::min(pixels[index  ]*255.f, 255.f), 0.f);
			unsigned char g = std::max(std::min(pixels[index+1]*255.f, 255.f), 0.f);
			unsigned char b = std::max(std::min(pixels[index+2]*255.f, 255.f), 0.f);
			fprintf(image, "%c%c%c", r,g,b);
		}
	}
	fclose(image);

	// Cleanup
	delete[] pixels;
}
