
#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/subsurface.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sensor.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/render/skdtree.h>
#include <mitsuba/render/scene.h>

#include <set>
#include <vector>

#include "grainmesh/sphere_pack.h"
#include "grainmesh/sphere_pack_data.h"

#include "grainmesh/grainmesh_kdtree.h"
#include "grainmesh/bvhinterface.h"

#include "grainmesh/tempinfo.h"


typedef unsigned char byte;

MTS_NAMESPACE_BEGIN

#define OUTPUT_CONDITION (control_timer->getMicroseconds()%10000 == 0)

#define TEST 0


#if !TEST


struct SphereIntersectionInfo{
    FastBVH::IntersectionInfo I;
    SpherePack::Sphere* sphere;
    int index;
    // mesh
};


class SpherePackDict{
public:
    SpherePackDict() {
        control_timer = new Timer();

        m_dictCount = sizeof(spherePacks) / sizeof(spherePacks[0]);
        m_sphereCount = sizeof(spherePacks[0]) / sizeof(struct SpherePack::Sphere);

        for(int i=0; i<m_dictCount; i++){
            BVHInterface* bvhInterface = new BVHInterface((struct SpherePack::Sphere*)spherePacks[i], m_sphereCount);
            m_spherePacks.push_back(bvhInterface);
        }
    }

    inline bool rayIntersect(int coord, const Ray& ray, SphereIntersectionInfo &info) {
        // which dict = f(x, y, z)
        int dictIndex = getDictType(coord);

        SpherePack::Sphere* sphere = (m_spherePacks.at(dictIndex))
                ->rayIntersect(ray.o, Point(ray.d.x, ray.d.y, ray.d.z), info.I);

        if(sphere){
            info.sphere = sphere;
            info.index = sphere - spherePacks[dictIndex]; // index in this pack
        }

        return NULL != sphere;
    }

    int getDictType(int coord) {
        srand(coord);
        return rand() % m_dictCount;
    }

    int getGrainType(int coord, int sphereIndex) {
        long res = 0;
        srand(coord);
        res += rand();
        srand(getDictType(coord));
        res += rand();
        srand(sphereIndex);
        res += rand();
        return res % m_grainCount;
    }

    Transform getGrainTransform(int coord, int sphereIndex) {
        int rotate;
        Transform tr;
        srand(coord);
        srand(rand() + getDictType(coord));
        srand(rand() + sphereIndex);

        // x
        rotate = rand() % 360;
        tr = tr.rotate(Vector(1,0,0), rotate);
        // y
        rotate = rand() % 360;
        tr = tr.rotate(Vector(0,1,0), rotate);
        // z
        rotate = rand() % 360;
        tr = tr.rotate(Vector(0,0,1), rotate);

        return tr;
    }

    void setGrainCount(int count) { m_grainCount = count; }

private:
    ref<Timer> control_timer;

    int m_dictCount;
    int m_sphereCount;
    int m_grainCount;
    std::vector<BVHInterface*> m_spherePacks;
};


class Grain {
public:
    grain::KDNode* m_kdtree;

    Grain(std::vector<grain::KDTriangle*>& triangels) {
        for(size_t i=0; i<triangels.size(); i++) {
            m_triangles.push_back(triangels.at(i));
        }

        // get a bounding box surrounding all the triangles
        AABB aabb = (triangels.at(0))->getAABB();
        for(size_t i=1; i < triangels.size(); i++) {
            aabb.expandBy((triangels.at(i))->getAABB());
        }

        // find the furthest vertex from sphere center
        Point midPoint = (aabb.min + aabb.max)/2.0;
        Float minD = 0.0, tmpD;

        for(size_t i=0; i<triangels.size(); i++) {
            for(int j=0; j<3; j++) {
                tmpD = distance(midPoint, ((grain::KDTriangle*)triangels.at(i))->m_points[j]);
                if(tmpD > minD)
                    minD = tmpD;
            }
        }

        // scale to a box inscribed inside a unit sphere
        Point translate = (aabb.min + aabb.max)/(-2.0);
        Float scale = 1.0 / minD;

        for (size_t i=0; i < triangels.size(); i++) {
            Point* p = ((grain::KDTriangle*)triangels.at(i))->m_points;
            for(size_t j=0; j<3; j++) {
                p[j] = (p[j] + translate) * scale;
            }
        }

        m_kdtree = new grain::KDNode(m_triangles);
        //m_aabb = m_kdtree->bbox;
    }

    void precomputeTSDF() {
        const int SAMPLE_NUM = 10;
        long hitCount = 0;

        Ray ray;
        ray.o = Point(-1.0, 0, 0);

        Float coord_theta, coord_phi, theta, phi;
        Vector coord_x, coord_y, coord_z;
        hitInfo hitinfo;

std::cout<<"[";

        // uniformly sample origin, on sphere
        for(int i=0; i<SAMPLE_NUM; i++) {
            coord_theta = (2.0*random()/RAND_MAX - 1.0) * M_PI; // -pi ~ pi
            coord_phi = (1.0*random()/RAND_MAX - 0.5) * M_PI; // -pi/2 ~ pi/2
            ray.o = (Point(cos(coord_phi)*cos(coord_theta), sin(coord_phi), -cos(coord_phi)*sin(coord_theta)));
            coord_x = Point(0,0,0) - ray.o;
            coord_z = normalize(cross(coord_x, Vector(0,1,0)));
            coord_y = normalize(cross(coord_z, coord_x));

            Matrix4x4 m(coord_x.x, coord_y.x, coord_z.x, 0,
                        coord_x.y, coord_y.y, coord_z.y, 0,
                        coord_x.z, coord_y.z, coord_z.z, 0,
                        0, 0, 0, 1);

            // uniformly sample direction, on hemi-sphere
            for(int j=0; j<SAMPLE_NUM; j++) {
                theta = (1.0*random()/RAND_MAX - 0.5) * M_PI; // -pi/2 ~ pi/2
                phi = (1.0*random()/RAND_MAX - 0.5) * M_PI; // -pi/2 ~ pi/2
                ray.d = normalize(Vector(cos(phi)*cos(theta), sin(phi), -cos(phi)*sin(theta)));
                Vector4 tmpd = m*Vector4(ray.d.x, ray.d.y, ray.d.z, 1);
                ray.d = normalize(Vector(tmpd.x, tmpd.y, tmpd.z));

                float tout=1.0/0.0, tmin=-0.;
                Float t;
                if(! m_kdtree->hit(ray, tout, tmin, &hitinfo))
                    continue;

                // hit
                t = (Float)tout;
                hitCount++;

                Normal n = static_cast<grain::KDTriangle*>(hitinfo.hitObject)->getNormal();
                Point hit = ray.o + t*ray.d;
                Vector outd = normalize(-2 * dot(ray.d, normalize(n)) * normalize(n) + ray.d);

                // cos(theta)
                Float cosTheta = dot(ray.d, outd);

                double A = outd.x*outd.x + outd.y*outd.y + outd.z*outd.z;
                double B = 2*(hit.x*outd.x + hit.y*outd.y + hit.z*outd.z);
                double C = hit.x*hit.x + hit.y*hit.y + hit.z*hit.z - 1.0;
                double nearT, farT;

                if (! solveQuadraticDouble(A, B, C, nearT, farT))
                    continue;

                Point out = hit + ((Float)farT) * outd;
                Point projection = ray.o + dot(ray.d, out-ray.o) * ray.d;

                // r, z
                Float r = distance(out, projection);
                Float z = distance(ray.o, projection);

                // phi
                Vector outdProjection = normalize(cross(cross(coord_x, outd), coord_x));
                double phi = acos((double)dot(outdProjection, Vector(0,1,0))) * 180 / M_PI;
                if(dot(coord_z, outdProjection) < 0)
                    phi = 360 - phi;


                //std::cout<<"("<<(float)cosTheta<<", "<<r<<", "<<z<<", "<<phi<<"), ";
            }
        }

        std::cout<<"]"<<std::endl;
        std::cout<<"hit probability: "<<hitCount * 1.0 / (SAMPLE_NUM*SAMPLE_NUM)<<std::endl;

    }

    inline bool rayIntersect(Transform tr, const Ray &_ray, hitInfo *info) {
        float tout = 1.0/0.0, tmin = 0.;
        Ray ray = tr.transformAffine(_ray);
        bool res = m_kdtree->hit(ray, tout, tmin, info);
        if(res) {
            info->hitPoint = tr.inverse().transformAffine(ray.o + ((Float)tout) * ray.d);
            return true;
        }
        return false;
    }

private:
    std::vector<grain::KDTriangle*> m_triangles;
    //AABB m_aabb;
};




class GrainMesh : public Shape {
private:
    Transform m_objectToWorld;
    Transform m_worldToObject;
    // TODO remove cylinder related
    Float m_radius, m_length, m_invSurfaceArea;
    int m_aggregateMeshCount, m_grainMeshCount;

    // voxel
    byte *m_voxels;
    int m_dimensionX, m_dimensionY, m_dimensionZ; // actually they are the same
    Point volumeOrigin; // origin of the voxels
    Vector volumeSpan; // whole volume size, cubical
    Vector gridSpan; // one voxel size

    // kd-tree
    grain::KDNode *m_mesh_kdtree;
    SpherePackDict *m_spherePackDict;


    std::vector<Grain*> m_grains;
    //std::vector<grain::KDTriangle*> m_grainTriangles;
    //grain::KDNode *m_grain_kdtree;

    // modes
    MODES_T m_mode;

    // debug
    ref<Timer> control_timer;
    bool simpledebug;


public:
    int read_binvox(std::string filespec,
                    bool first=true /*Is this first time to read? if yes, need to allocate memory*/)
    {
        fs::ifstream *input = new fs::ifstream(filespec.c_str(), std::ios::in | std::ios::binary);
        if (input->bad() || input->fail())
            Log(EError, "Binvox file '%s' not found!", filespec.c_str());

        //
        // read header
        //
        std::string line;
        *input >> line;  // #binvox
        if (line.compare("#binvox") != 0) {
            Log(EError, "Error: first line reads [%s] instead of [#binvox]", line.c_str());
            delete input;
            return 0;
        }
        int version;
        *input >> version;
        Log(EInfo, "reading binvox version %d", version);

        int depth, height, width;
        depth = -1;
        int done = 0;
        while(input->good() && !done) {
            *input >> line;
            if (line.compare("data") == 0) done = 1;
            else if (line.compare("dim") == 0) {
                *input >> depth >> height >> width;
                if(first) {
                    m_dimensionX = width;
                    m_dimensionY = height;
                    m_dimensionZ = depth;
                }
            }
            else if (line.compare("translate") == 0) {
                float tx, ty, tz;
                *input >> tx >> ty >> tz;
                if(first) {
                    volumeOrigin.x = tx;
                    volumeOrigin.y = ty;
                    volumeOrigin.z = tz;
                }
            }
            else if (line.compare("scale") == 0) {
                float size;
                *input >> size;
                if(first) {
                    volumeSpan.x = volumeSpan.y = volumeSpan.z = size;
                }
            }
            else {
                Log(EWarn, "  unrecognized keyword [%s], skipping", line.c_str());
                char c;
                do {  // skip until end of line
                    c = input->get();
                } while(input->good() && (c != '\n'));

            }
        }
        if (!done) {
            Log(EError, "  error reading header");
            return 0;
        }
        if (depth == -1) {
            Log(EError, "  missing dimensions in header");
            return 0;
        }

        //
        // allocate and init
        //
        int size = width * height * depth;
        if(first) {
            m_voxels = new byte[size];
            if (!m_voxels) {
                Log(EError, "  error allocating memory");
                return 0;
            }
            for(int i=0; i<size; i++)
                m_voxels[i] = 0;
        }

        //
        // read voxel data
        //
        byte value;
        byte count;
        int index = 0;
        int end_index = 0;
        int nr_voxels = 0;

        input->unsetf(std::ios::skipws);  // need to read every byte now (!)
        *input >> value;  // read the linefeed char

        while((end_index < size) && input->good()) {
            *input >> value >> count;

            if (input->good()) {
                end_index = index + count;
                if (end_index > size) return 0;
                for(int i=index; i < end_index; i++)
                    m_voxels[i] = m_voxels[i] + value;

                if (value) nr_voxels += count;
                index = end_index;
            }  // if file still ok

        }  // while

        input->close();
        Log(EInfo, "  read %d voxels", nr_voxels);

        return 1;

    }

    void extractKDTrianglesFromMeshes(
            std::vector<TriMesh*> &meshes,
            std::vector<grain::KDTriangle*> &KDTriangles)
    {
        for(size_t i=0; i<meshes.size(); i++)
        {
            int count = meshes.at(i)->getTriangleCount();
            const Triangle* triangles = meshes.at(i)->getTriangles();
            Point* positions = meshes.at(i)->getVertexPositions();
            Normal* normals;
            if(meshes.at(i)->hasVertexNormals())
                normals = meshes.at(i)->getVertexNormals();
            Point2* uvs;
            if(meshes.at(i)->hasVertexTexcoords())
                uvs = meshes.at(i)->getVertexTexcoords();

            for(size_t j=0; j<count; j++) {
                grain::KDTriangle *t = new grain::KDTriangle(
                            positions[triangles[j].idx[0]],
                            positions[triangles[j].idx[1]],
                            positions[triangles[j].idx[2]]);

                if(meshes.at(i)->hasVertexNormals()) {
                    t->addNormal(normals[triangles[j].idx[0]],
                                 normals[triangles[j].idx[1]],
                                 normals[triangles[j].idx[2]]);
                }

                if(meshes.at(i)->hasVertexTexcoords()) {
                    t->addUV(uvs[triangles[j].idx[0]],
                             uvs[triangles[j].idx[1]],
                             uvs[triangles[j].idx[2]]);
                }

                KDTriangles.push_back(t);
            }
        }
    }

    GrainMesh(const Properties &props) : Shape(props) {

        control_timer = new Timer();
        m_radius = 1.0, m_length = 2.0;

        ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver()->clone();


        // stub
        m_aggregateMeshCount = 1;
        fs::path path = fileResolver->resolve(props.getString("filename"));
        m_name = path.stem().string();

        /* By default, any existing normals will be used for
           rendering. If no normals are found, Mitsuba will
           automatically generate smooth vertex normals.
           Setting the 'faceNormals' parameter instead forces
           the use of face normals, which will result in a faceted
           appearance.
        */
        m_faceNormals = props.getBoolean("faceNormals", false);

        /* Causes all normals to be flipped */
        m_flipNormals = props.getBoolean("flipNormals", false);

        /* Collapse all contained shapes / groups into a single object? */
        m_collapse = props.getBoolean("collapse", false);

        /* Causes all texture coordinates to be vertically flipped */
        bool flipTexCoords = props.getBoolean("flipTexCoords", true);

        /// When the file contains multiple meshes, this index specifies which one to load
        int shapeIndex = props.getInteger("shapeIndex", -1);

        /* Object-space -> World-space transformation */
        Transform objectToWorld = props.getTransform("toWorld", Transform());
        m_objectToWorld = objectToWorld;
        m_worldToObject = m_objectToWorld.inverse();

        /* Import materials from a MTL file, if any? */
        bool loadMaterials = props.getBoolean("loadMaterials", true);

        /* Load the geometry */
        Log(EInfo, "Loading geometry from \"%s\" ..", path.filename().string().c_str());
        fs::ifstream is(path);
        if (is.bad() || is.fail())
            Log(EError, "Wavefront OBJ file '%s' not found!", path.string().c_str());

        fileResolver->prependPath(fs::absolute(path).parent_path());

        ref<Timer> timer = new Timer();

        loadWavefrontOBJ(m_meshes, m_materialAssignment,
                         m_name, is, fileResolver, m_objectToWorld, m_collapse, flipTexCoords, loadMaterials, shapeIndex);

        /* smooth angle */
        if (props.hasProperty("maxSmoothAngle")) {
            if (m_faceNormals)
                Log(EError, "The properties 'maxSmoothAngle' and 'faceNormals' "
                "can't be specified at the same time!");
            Float maxSmoothAngle = props.getFloat("maxSmoothAngle");
            for (size_t i=0; i<m_meshes.size(); ++i)
                m_meshes[i]->rebuildTopology(maxSmoothAngle);
        }

        Log(EInfo, "Done with \"%s\" (took %i ms)", path.filename().string().c_str(), timer->getMilliseconds());
        timer->reset();


        // Create the mesh KD-Tree for intersection
        std::vector<grain::KDTriangle*> allKDTriangles;
        extractKDTrianglesFromMeshes(m_meshes, allKDTriangles);

        Log(EInfo, "KD-tree: %d triangles in total", allKDTriangles.size());
        m_mesh_kdtree = new grain::KDNode(allKDTriangles);

        Log(EInfo, "Done build the kd-tree (took %i ms)", timer->getMilliseconds());
        timer->reset();

        /* Load the binvox file */
        Log(EInfo, "Loading binvox from \"%s\" ..", path.filename().string().c_str());

        read_binvox(path.string()+".binvox");
        read_binvox(path.string()+".hollow.binvox", false);

        m_volume_aabb.min = volumeOrigin;
        m_volume_aabb.max = volumeOrigin + volumeSpan;

        Log(EInfo, "Loaded binvox with \"%s\" (took %i ms)", path.filename().string().c_str(), timer->getMilliseconds());
        timer->reset();

        /* Load the instantiation dictionary */
        m_spherePackDict = new SpherePackDict();

        Log(EInfo, "Loaded instantiation dictionary (took %i ms)", timer->getMilliseconds());
        timer->reset();


        /* Load grain meshes */
        m_grainMeshCount = props.getInteger("graintype", 0);
        if(m_grainMeshCount < 1)
            Log(EError, "No grain mesh!");

        m_spherePackDict->setGrainCount(m_grainMeshCount);

        for (int i=0; i<m_grainMeshCount; i++) {
            std::stringstream ss;
            ss<<"grain"<<i;

            fs::path path = fileResolver->resolve(props.getString(ss.str()));

            Log(EInfo, "Loading grain %d from \"%s\" ..", i, path.filename().string().c_str());
            fs::ifstream is(path);
            if (is.bad() || is.fail())
                Log(EError, "Wavefront OBJ file '%s' not found!", path.string().c_str());

            fileResolver->prependPath(fs::absolute(path).parent_path());

            std::vector<TriMesh*> meshes;
            std::vector<std::string> materials;
            loadWavefrontOBJ(meshes, materials, ss.str(), is, fileResolver, m_objectToWorld);


            // create kd-tree for this grain
            std::vector<grain::KDTriangle*> grainTriangles;
            extractKDTrianglesFromMeshes(meshes, grainTriangles);
            m_grains.push_back(new Grain(grainTriangles));
        }

        Log(EInfo, "Loaded %d grain mesh(es) (took %i ms)", m_grainMeshCount, timer->getMilliseconds());
        timer->reset();

        /* precompute TSDF */
        ((Grain*)m_grains.at(0))->precomputeTSDF();
        //for(size_t i=0; i<m_grains.size(); i++) {
          //  ((Grain*)m_grains.at(i))->precomputeTSDF();
        //}

        Log(EInfo, "Precomputed TSDFs for %d grain mesh(es) (took %i ms)",m_grainMeshCount, timer->getMilliseconds());
        timer->reset();


        // TODO mitsuba::ShapeKDTree kdt;


        //==========debug=tweak==================
        int rendermode = props.getInteger("rendermode", 0);
        switch (rendermode) {
        case 0:
            m_mode = EPT;
            break;
        case 1:
            m_mode = VPT;
            break;
        case 2:
            m_mode = DA;
            break;
        default:
            m_mode = DEBUG;
            break;
        }


        simpledebug = props.getBoolean("simpledebug", false);
        if(simpledebug) {
            for(int i=0; i<8; i++)
                m_voxels[i] = 0;

            m_voxels[0] = 1;
            m_voxels[1] = 1;
            m_voxels[2] = 1;
            m_voxels[4] = 1;
            m_dimensionX = 2;
            m_dimensionY = 2;
            m_dimensionZ = 2;
        }

    }

    void configure() {
        Shape::configure();

        m_aabb.reset();
        for (size_t i=0; i<m_meshes.size(); ++i) {
            m_meshes[i]->configure();
            m_aabb.expandBy(m_meshes[i]->getAABB());
        }

        // volumeSpan = m_aabb.max - m_aabb.min;
        gridSpan.x = volumeSpan.x / m_dimensionX;
        gridSpan.y = volumeSpan.y / m_dimensionY;
        gridSpan.z = volumeSpan.z / m_dimensionZ;

        Log(EInfo, "Volume span: %f %f %f", volumeSpan.x, volumeSpan.y, volumeSpan.z);
        Log(EInfo, "Grid span: %f %f %f", gridSpan.x, gridSpan.y, gridSpan.z);
        Log(EInfo, "mesh aabb:\n%s", m_aabb.toString().c_str());
        Log(EInfo, "volume aabb:\n%s", m_volume_aabb.toString().c_str());




#if 0
        AABB aabb = m_volume_aabb;

        size_t size = m_dimensionX*m_dimensionY*m_dimensionZ;
        byte *inval = new byte[size];
        for(size_t i=0; i<size; i++)
            inval[i] = 0;

        int index;
        for(int x=0; x<m_dimensionX; x++)
            for(int y=0; y<m_dimensionY; y++)
                for(int z=0; z<m_dimensionZ; z++) {

                    Point midp(aabb.min.x+gridSpan.x*(x+0.5), aabb.min.y+gridSpan.y*(y+0.5), aabb.min.z+gridSpan.z*(z+0.5));
                    bool res = m_mesh_kdtree->isPointInside(midp);

                    if(res) {
                        index = x * m_dimensionY * m_dimensionZ + z * m_dimensionY + y;
                        inval[index] = inval[index] + 1;
                    }
                }

        int err = 0;
        int expp = 0;
        for(size_t i=0; i<size; i++) {
            if((m_voxels[i]==0 && inval[i]==0) || (m_voxels[i]==1 && inval[i]==1))
                ;
            else
                err++;

            if(m_voxels[i]==2)
                expp++;
        }


        std::cout<<"error count: "<<err<<", "<<err*1.0/(size-expp)<<"%"<<std::endl;

#endif




        // TODO determine EPT, VPT, DA
#if 0
        const Sensor * sen = getSensor();
        if(sen) {
            Vector4 ori= ((const AnimatedTransform*)sen->getWorldTransform())->getTransform().getMatrix().col(3);
            Point cameraOrigin = Point(ori.x, ori.y, ori.z);

            if(distance(cameraOrigin, (m_volume_aabb.min+m_volume_aabb.max)/2) < 4000) {
                Log(EWarn, "Using EPT.");
                m_mode = EPT;
            }
            else {
                Log(EWarn, "Using VPT.");
                m_mode = VPT;
            }

            {
                m_mode = DEBUG;
                Log(EWarn, "debug, switch to DEBUG");
            }
        }
#endif




#if 0

        {
            std::vector<grain::KDTriangle*> tris;

            tris.push_back(new grain::KDTriangle(
                               Point(-1,-1,1),
                               Point(1,-1,1),
                               Point(0,1,1)
                               ));

            grain::KDNode* kdtree = new grain::KDNode(tris);

            Point o(0,0,2);
            Vector d(0,0,-1);
            Ray ray = Ray(o, d, (Float)0);

            float t=1.0/0.0, tmin=0;
            bool res = kdtree->hit(ray, t, tmin);
            if(res)
                Log(EWarn, "tttttttttttying: %f", t);
            if(res && t>tmin)
                Log(EWarn, "tttttttttttying: hit");
            else
                Log(EWarn, "tttttttttttying: not hit");

            //exit(0);

        }
#endif

    }

    GrainMesh(Stream *stream, InstanceManager *manager) : Shape(stream, manager) {
        m_aabb = AABB(stream);
        uint32_t meshCount = stream->readUInt();
        m_meshes.resize(meshCount);

        for (uint32_t i=0; i<meshCount; ++i) {
            m_meshes[i] = static_cast<TriMesh *>(manager->getInstance(stream));
            m_meshes[i]->incRef();
        }
    }

    bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {

        // Transform into the local coordinate system and normalize
        Ray ray;
        m_worldToObject(_ray, ray);

        // EPT
        if(EPT == m_mode){

            Point rayOrigin = ray.o;

            Vector dirinv;
            dirinv.x = 1.0f / ray.d.x;
            dirinv.y = 1.0f / ray.d.y;
            dirinv.z = 1.0f / ray.d.z;

            AABB aabb = m_volume_aabb;

            // if ray origin is outside, find the first intersection
            if(ray.o.x < aabb.min.x || ray.o.x >= aabb.max.x
            || ray.o.y < aabb.min.y || ray.o.y >= aabb.max.y
            || ray.o.z < aabb.min.z || ray.o.z >= aabb.max.z)
            {
                Float t1 = (aabb.min.x - ray.o.x) * dirinv.x;
                Float t2 = (aabb.max.x - ray.o.x) * dirinv.x;
                Float t3 = (aabb.min.y - ray.o.y) * dirinv.y;
                Float t4 = (aabb.max.y - ray.o.y) * dirinv.y;
                Float t5 = (aabb.min.z - ray.o.z) * dirinv.z;
                Float t6 = (aabb.max.z - ray.o.z) * dirinv.z;

                Float tmin = fmax(fmax(fmin(t1, t2), fmin(t3, t4)), fmin(t5, t6));
                Float tmax = fmin(fmin(fmax(t1, t2), fmax(t3, t4)), fmax(t5, t6));

                // if tmax < 0, ray (line) is intersecting AABB, but whole AABB is behing us
                if (tmax < 0) {
                    return false;
                }
                // if tmin > tmax, ray doesn't intersect AABB
                if (tmin > tmax) {
                    return false;
                }

                // t = tmin; return true;

                //Float tPush = (gridSpan.x * 0.01) / fmax(fmax(fabs(ray.d.x), fabs(ray.d.y)), fabs(ray.d.z)); // to push it inside the cube
                //tPush = 0.0;
                ray.o = ray.o + ray.d * (tmin);
            }

            // march in voxels to find an intersection
            int x, y, z;
            int stepX=0, stepY=0, stepZ=0;
            double tXMax=1.0/0.0, tYMax=1.0/0.0, tZMax=1.0/0.0; // inf
            double tXDelta=0, tYDelta=0, tZDelta=0;

            x = (int)((ray.o.x - aabb.min.x) / gridSpan.x);
            if(x == m_dimensionX) x-=1; // on the max border exactly
            y = (int)((ray.o.y - aabb.min.y) / gridSpan.y);
            if(y == m_dimensionY) y-=1;
            z = (int)((ray.o.z - aabb.min.z) / gridSpan.z);
            if(z == m_dimensionZ) z-=1;

            //if(x < 0 || x >= m_dimensionX || y < 0 || y >= m_dimensionY || z < 0 || z >= m_dimensionZ)
                //return false;

            const double TINY_FLOAT = 1e-20;
            stepX = ray.d.x > TINY_FLOAT ? 1 : (ray.d.x < -TINY_FLOAT ? -1 : 0);
            stepY = ray.d.y > TINY_FLOAT ? 1 : (ray.d.y < -TINY_FLOAT ? -1 : 0);
            stepZ = ray.d.z > TINY_FLOAT ? 1 : (ray.d.z < -TINY_FLOAT ? -1 : 0);
            if(stepX!=0) {
                tXMax = (gridSpan.x * (stepX>0 ? (x+1) : x) + aabb.min.x - ray.o.x) * dirinv.x;
                tXDelta = gridSpan.x * dirinv.x * stepX;
            }
            if(stepY!=0) {
                tYMax = (gridSpan.y * (stepY>0 ? (y+1) : y) + aabb.min.y - ray.o.y) * dirinv.y;
                tYDelta = gridSpan.y * dirinv.y * stepY;
            }
            if(stepZ!=0) {
                tZMax = (gridSpan.z * (stepZ>0 ? (z+1) : z) + aabb.min.z - ray.o.z) * dirinv.z;
                tZDelta = gridSpan.z * dirinv.z * stepZ;
            }


            SphereIntersectionInfo info;
            Vector inc(0,0,0);
            while(1) {
                bool test_sphere = true;

                // test voxel. voxel empty? full? half-empty?
                int voxelIndex = x * m_dimensionY * m_dimensionZ + z * m_dimensionY + y;
                int occupy = (int)m_voxels[voxelIndex];
                if(occupy == 0)
                    test_sphere = false;

                while(test_sphere) {
                    // test sphere
                    if(! m_spherePackDict->rayIntersect(
                            voxelIndex,
                            // transform into the unit voxel coordinate
                            Ray(Point((ray.o.x - aabb.min.x) / gridSpan.x - x,
                                      (ray.o.y - aabb.min.y) / gridSpan.y - y,
                                      (ray.o.z - aabb.min.z) / gridSpan.z - z),
                                ray.d, ray.time),
                            info))
                        break;


                    const FastBVH::Sphere* sphere = static_cast<const FastBVH::Sphere*>(info.I.object);
                    Point center;
                    center.x = (sphere->center.x + x) * gridSpan.x + aabb.min.x;
                    center.y = (sphere->center.y + y) * gridSpan.y + aabb.min.y;
                    center.z = (sphere->center.z + z) * gridSpan.z + aabb.min.z;

                    bool test_grain = true;
                    // half-empty
                    if(occupy == 2) {
                        if(m_mesh_kdtree->isPointInside(center)) {
                            ; //test_grain = true;
                        }
                        else
                            test_grain = false;
                    }

                    // test grain mesh
                    if(test_grain) {
                        hitInfo hitinfo;
                        Transform transform = m_spherePackDict->getGrainTransform(voxelIndex, info.index);
                        bool res = ((Grain*)m_grains.at(m_spherePackDict->getGrainType(voxelIndex, info.index)))->rayIntersect(
                                    transform,
                                    // transform into the unit sphere
                                    Ray(Point(ray.o.x-center.x, ray.o.y-center.y, ray.o.z-center.z) / (sphere->r * gridSpan.x), ray.d, ray.time),
                                    &hitinfo
                                    );

                        if(res) {
                            Point hitPoint;
                            hitPoint.x = hitinfo.hitPoint.x * gridSpan.x * sphere->r + center.x;
                            hitPoint.y = hitinfo.hitPoint.y * gridSpan.y * sphere->r + center.y;
                            hitPoint.z = hitinfo.hitPoint.z * gridSpan.z * sphere->r + center.z;

                            // actually, they are almost the same
                            t = fmin(fmin((hitPoint.x - rayOrigin.x) / ray.d.x, (hitPoint.y - rayOrigin.y) / ray.d.y), (hitPoint.z - rayOrigin.z) / ray.d.z);

                            if(temp){
                                // init the temp space, MTS_KD_INTERSECTION_TEMP-4 bytes
                                new(temp) GrainIntersectionInfo;
                                GrainIntersectionInfo* tempinfo = static_cast<GrainIntersectionInfo*>(temp);
                                tempinfo->valid = true;
                                tempinfo->rendermode = EPT;
                                tempinfo->object = hitinfo.hitObject;
                                tempinfo->tr = transform;
                            }

                            return true;
                        }
                    }

                    // skip this sphere. sphere not in aggregate mesh or grain not hit.
                    Point hitPoint;
                    hitPoint.x = (info.I.hit.x + x) * gridSpan.x + aabb.min.x;
                    hitPoint.y = (info.I.hit.y + y) * gridSpan.y + aabb.min.y;
                    hitPoint.z = (info.I.hit.z + z) * gridSpan.z + aabb.min.z;

                    ray.o = hitPoint +
                            (2 * (dot(ray.d,
                                      (Vector(center.x, center.y, center.z)
                                       - Vector(hitPoint.x, hitPoint.y, hitPoint.z)))))
                            * ray.d;

                } // end of sphere test while


                // perform intersection above, do not miss the first voxel
                // =====================================

                // find_next_voxel:

                if(tXMax <= tYMax && tXMax <= tZMax) { //x
                    x += stepX;
                    if(x < 0 || x >= m_dimensionX)
                        return false;
                    inc.x = 1.0;
                    //tXMax += tXDelta;
                } else {
                    inc.x = 0.0;
                }

                if(tYMax <= tXMax && tYMax <= tZMax) { // y
                    y += stepY;
                    if(y < 0 || y >= m_dimensionY)
                        return false;
                    inc.y = 1.0;
                    //tYMax += tYDelta;
                } else {
                    inc.y = 0.0;
                }
                if(tZMax <= tXMax && tZMax <= tYMax) { //z
                    z += stepZ;
                    if(z < 0 || z >= m_dimensionZ)
                        return false;
                    inc.z = 1.0;
                    //tZMax += tZDelta;
                } else {
                    inc.z = 0.0;
                }

                tXMax += tXDelta * inc.x;
                tYMax += tYDelta * inc.y;
                tZMax += tZDelta * inc.z;
            }

            return false;
        }
        // VPT
        else if(VPT == m_mode) {

            if(temp) {
                new(temp) GrainIntersectionInfo;
                GrainIntersectionInfo* tempinfo = static_cast<GrainIntersectionInfo*>(temp);
                tempinfo->valid = false;
                tempinfo->rendermode = VPT;
            }

            Point center = ray.o;
            if(m_mesh_kdtree->isPointInside(center))
            {
                double theta_t = 1.0 / (0.87 + 1.0 * 0.18);
                double v = 1.0*random()/RAND_MAX;
                t = log(1-v) / (-theta_t) * gridSpan.x / 2;
                return true;
            }
            else {
                float tout = 1.0/0.0, tmin = 0.;
                bool res = m_mesh_kdtree->hit(ray, tout, tmin);
                if(res) {
                    t = (Float)tout;
                    return true;
                }
                return false;
            }
        }
        // DA
        else if(DA == m_mode) {
            ;
        }
        // debug
        else {
            hitInfo hitinfo;
            bool res = ((Grain*)m_grains.at(0))->rayIntersect( // TODO, grain index
                        Transform(),
                        // transform into the unit sphere
                        Ray(Point(ray.o.x, ray.o.y - gridSpan.y, ray.o.z) / (500.0), ray.d, ray.time),
                        &hitinfo
                        );

            if(res) {
                t = 1.0;

                if(temp) {
                    new(temp) GrainIntersectionInfo;
                    GrainIntersectionInfo* tempinfo = (GrainIntersectionInfo*)temp;
                    tempinfo->valid = true;
                    tempinfo->rendermode = DEBUG;
                    tempinfo->object = hitinfo.hitObject;
                }

                return true;
                }
            return false;
        }

        return false;
    }

    bool rayIntersect(const Ray &_ray, Float mint, Float maxt) const {
        Ray ray;
        m_worldToObject(_ray, ray);

        if(EPT == m_mode) {
            Float t;
            return rayIntersect(_ray, mint, maxt, t, NULL);
        }
        else if(VPT == m_mode) {
            ;
        }
        else if(DA == m_mode) {
            ;
        }
        else {
            ;
        }

        /*
        float tout;
        return m_mesh_kdtree->hit(ray, tout, mint);
        */

        return false;        
    }

    void fillIntersectionRecord(const Ray &ray,
            const void *temp, Intersection &its) const {
        its.p = ray(its.t);
        Point local = m_worldToObject(its.p);

        const GrainIntersectionInfo* tempinfo = static_cast<const GrainIntersectionInfo*>(temp);

        if(tempinfo->valid) {
            grain::KDTriangle* kdtriangle = static_cast<grain::KDTriangle*>(tempinfo->object);

            Transform tr = tempinfo->tr.inverse();

            its.uv = kdtriangle->getUV();
            its.shape = this;
            TangentSpace tangentSpace = kdtriangle->getUVTangent(); // TODO transform
            its.dpdu = m_objectToWorld(tr(tangentSpace.dpdu));
            its.dpdv = m_objectToWorld(tr(tangentSpace.dpdv));
            its.geoFrame.s = normalize(its.dpdu);
            its.geoFrame.t = normalize(its.dpdv);
            its.geoFrame.n = Normal(cross(its.geoFrame.s, its.geoFrame.t));

            /* Mitigate roundoff error issues by a normal shift of the computed intersection point */
            //its.p += its.geoFrame.n * (m_radius - std::sqrt(local.x*local.x+local.y*local.y));

            if (m_flipNormals)
                its.geoFrame.n *= -1;
            its.shFrame.n = its.geoFrame.n;
            its.hasUVPartials = false;
            its.instance = NULL;
            its.time = ray.time;


            // TODO eliminate temp
            its.temp = (void*)temp;
        }
        else {
            Float phi = std::atan2(local.y, local.x);
            if (phi < 0)
                phi += 2*M_PI;
            its.uv.x = phi / (2*M_PI);
            its.uv.y = local.z / m_length;

            Vector dpdu = Vector(-local.y, local.x, 0) * (2*M_PI);
            Vector dpdv = Vector(0, 0, m_length);
            its.shape = this;
            its.dpdu = m_objectToWorld(dpdu);
            its.dpdv = m_objectToWorld(dpdv);
            its.geoFrame.s = normalize(its.dpdu);
            its.geoFrame.t = normalize(its.dpdv);
            its.geoFrame.n = Normal(cross(its.geoFrame.s, its.geoFrame.t));

            /* Migitate roundoff error issues by a normal shift of the computed intersection point */
            its.p += its.geoFrame.n * (m_radius - std::sqrt(local.x*local.x+local.y*local.y));

            if (m_flipNormals)
                its.geoFrame.n *= -1;
            its.shFrame.n = its.geoFrame.n;
            its.hasUVPartials = false;
            its.instance = NULL;
            its.time = ray.time;
        }


        /*
        Float vvv[SPECTRUM_SAMPLES];
        vvv[0] = 0.5;
        vvv[1] = 0.2;
        vvv[2] = 0.3;
        Spectrum value(vvv);
        its.color = value;
        */
    }

    struct OBJTriangle {
        int p[3];
        int n[3];
        int uv[3];

        inline OBJTriangle() {
            memset(this, 0, sizeof(OBJTriangle));
        }
    };

    struct Vertex {
        Point p;
        Normal n;
        Point2 uv;
    };

    bool fetch_line(std::istream &is, std::string &line) {
        /// Fetch a line from the stream, while handling line breaks with backslashes
        if (!std::getline(is, line))
            return false;
        if (line == "")
            return true;
        int lastCharacter = (int) line.size() - 1;
        while (lastCharacter >= 0 &&
                (line[lastCharacter] == '\r' ||
                    line[lastCharacter] == '\n' ||
                    line[lastCharacter] == '\t' ||
                    line[lastCharacter] == ' '))
            lastCharacter--;

        if (lastCharacter >= 0 && line[lastCharacter] == '\\') {
            std::string nextLine;
            fetch_line(is, nextLine);
            line = line.substr(0, lastCharacter) + nextLine;
        } else {
            line.resize(lastCharacter+1);
        }
        return true;
    }

    void parse(OBJTriangle &t, int i, const std::string &str) {
        std::vector<std::string> tokens = tokenize(str, "/");
        if (tokens.size() == 1) {
            t.p[i] = atoi(tokens[0].c_str());
        } else if (tokens.size() == 2) {
            if (str.find("//") == std::string::npos) {
                t.p[i]  = atoi(tokens[0].c_str());
                t.uv[i] = atoi(tokens[1].c_str());
            } else {
                t.p[i] = atoi(tokens[0].c_str());
                t.n[i] = atoi(tokens[1].c_str());
            }
        } else if (tokens.size() == 3) {
            t.p[i] = atoi(tokens[0].c_str());
            t.uv[i] = atoi(tokens[1].c_str());
            t.n[i] = atoi(tokens[2].c_str());
        } else {
            Log(EError, "Invalid OBJ face format!");
        }
    }

    Texture *loadTexture(const FileResolver *fileResolver,
            std::map<std::string, Texture *> &cache,
            const fs::path &mtlPath, std::string filename,
            bool noGamma = false) {
        /* Prevent Linux/OSX fs::path handling issues for DAE files created on Windows */
        for (size_t i=0; i<filename.length(); ++i) {
            if (filename[i] == '\\')
                filename[i] = '/';
        }

        if (cache.find(filename) != cache.end())
            return cache[filename];

        fs::path path = fileResolver->resolve(filename);
        if (!fs::exists(path)) {
            path = fileResolver->resolve(fs::path(filename).filename());
            if (!fs::exists(path)) {
                Log(EWarn, "Unable to find texture \"%s\" referenced from \"%s\"!",
                    path.string().c_str(), mtlPath.string().c_str());
                return new ConstantSpectrumTexture(Spectrum(0.0f));
            }
        }
        Properties props("bitmap");
        props.setString("filename", path.string());
        if (noGamma)
            props.setFloat("gamma", 1.0f);
        ref<Texture> texture = static_cast<Texture *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Texture), props));
        texture->configure();
        texture->incRef();
        cache[filename] = texture;
        return texture;
    }

    void addMaterial(const std::string &name, Texture *diffuse, Texture *specular,
            Texture *exponent, Texture *bump, Texture *mask, int model) {
        ref<BSDF> bsdf;
        Properties props;

        if (model == 2 && (specular->getMaximum().isZero() || exponent->getMaximum().isZero()))
            model = 1;

        if (model == 2) {
            props.setPluginName("phong");

            bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(BSDF), props));
            bsdf->addChild("diffuseReflectance", diffuse);
            bsdf->addChild("specularReflectance", specular);
            bsdf->addChild("exponent", exponent);
        } else if (model == 4 || model == 6 || model == 7 || model == 9) {
            props.setPluginName("dielectric");
            bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(BSDF), props));
        } else if (model == 5 || model == 8) {
            props.setPluginName("conductor");
            props.setString("material", "Al");
            bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(BSDF), props));
        } else {
            props.setPluginName("diffuse");
            bsdf = static_cast<BSDF *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(BSDF), props));
            bsdf->addChild("reflectance", diffuse);
        }

        bsdf->configure();

        if (bump) {
            props = Properties("bumpmap");
            ref<BSDF> bumpBSDF = static_cast<BSDF *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(BSDF), props));

            bumpBSDF->addChild(bump);
            bumpBSDF->addChild(bsdf);
            bumpBSDF->configure();
            bsdf = bumpBSDF;
        }

        if (mask) {
            props = Properties("mask");
            ref<BSDF> maskedBSDF = static_cast<BSDF *> (PluginManager::getInstance()->
                createObject(MTS_CLASS(BSDF), props));

            maskedBSDF->addChild("opacity", mask);
            maskedBSDF->addChild(bsdf);
            maskedBSDF->configure();
            bsdf = maskedBSDF;
        }

        bsdf->setID(name);
        addChild(name, bsdf, false);
    }

    void loadMaterialLibrary(const FileResolver *fileResolver, const fs::path &mtlPath) {
        if (!fs::exists(mtlPath)) {
            Log(EWarn, "Could not find referenced material library '%s'",
                mtlPath.string().c_str());
            return;
        }

        Log(EInfo, "Loading OBJ materials from \"%s\" ..", mtlPath.filename().string().c_str());
        fs::ifstream is(mtlPath);
        if (is.bad() || is.fail())
            Log(EError, "Unexpected I/O error while accessing material file '%s'!",
                mtlPath.string().c_str());
        std::string buf, line;
        std::string mtlName;
        ref<Texture> specular, diffuse, exponent, bump, mask;
        int illum = 0;
        specular = new ConstantSpectrumTexture(Spectrum(0.0f));
        diffuse = new ConstantSpectrumTexture(Spectrum(0.0f));
        exponent = new ConstantFloatTexture(0.0f);
        std::map<std::string, Texture *> cache;

        while (is.good() && !is.eof() && fetch_line(is, line)) {
            std::istringstream iss(line);
            if (!(iss >> buf))
                continue;

            if (buf == "newmtl") {
                if (mtlName != "")
                    addMaterial(mtlName, diffuse, specular, exponent, bump, mask, illum);

                mtlName = trim(line.substr(6, line.length()-6));

                specular = new ConstantSpectrumTexture(Spectrum(0.0f));
                diffuse = new ConstantSpectrumTexture(Spectrum(0.0f));
                exponent = new ConstantFloatTexture(0.0f);
                mask = NULL;
                bump = NULL;
                illum = 0;
            } else if (buf == "Kd") {
                Float r, g, b;
                iss >> r >> g >> b;
                Spectrum value;
                value.fromSRGB(r, g, b);
                diffuse = new ConstantSpectrumTexture(value);
            } else if (buf == "map_Kd") {
                std::string filename;
                iss >> filename;
                diffuse = loadTexture(fileResolver, cache, mtlPath, filename);
            } else if (buf == "Ks") {
                Float r, g, b;
                iss >> r >> g >> b;
                Spectrum value;
                value.fromSRGB(r, g, b);
                specular = new ConstantSpectrumTexture(value);
            } else if (buf == "map_Ks") {
                std::string filename;
                iss >> filename;
                specular = loadTexture(fileResolver, cache, mtlPath, filename);
            } else if (buf == "bump") {
                std::string filename;
                iss >> filename;
                bump = loadTexture(fileResolver, cache, mtlPath, filename, true);
            } else if (buf == "map_d") {
                std::string filename;
                iss >> filename;
                mask = loadTexture(fileResolver, cache, mtlPath, filename);
            } else if (buf == "d" /* || buf == "Tr" */) {
                Float value;
                iss >> value;
                if (value == 1)
                    mask = NULL;
                else
                    mask = new ConstantFloatTexture(value);
            } else if (buf == "Ns") {
                Float value;
                iss >> value;
                exponent = new ConstantFloatTexture(value);
            } else if (buf == "illum") {
                iss >> illum;
            } else {
                /* Ignore */
            }
        }

        addMaterial(mtlName, diffuse, specular, exponent, bump, mask, illum);

        for (std::map<std::string, Texture *>::iterator it = cache.begin();
                it != cache.end(); ++it)
            it->second->decRef();
    }

    /// For using vertices as keys in an associative structure
    struct vertex_key_order : public
        std::binary_function<Vertex, Vertex, bool> {
    public:
        bool operator()(const Vertex &v1, const Vertex &v2) const {
            if (v1.p.x < v2.p.x) return true;
            else if (v1.p.x > v2.p.x) return false;
            if (v1.p.y < v2.p.y) return true;
            else if (v1.p.y > v2.p.y) return false;
            if (v1.p.z < v2.p.z) return true;
            else if (v1.p.z > v2.p.z) return false;
            if (v1.n.x < v2.n.x) return true;
            else if (v1.n.x > v2.n.x) return false;
            if (v1.n.y < v2.n.y) return true;
            else if (v1.n.y > v2.n.y) return false;
            if (v1.n.z < v2.n.z) return true;
            else if (v1.n.z > v2.n.z) return false;
            if (v1.uv.x < v2.uv.x) return true;
            else if (v1.uv.x > v2.uv.x) return false;
            if (v1.uv.y < v2.uv.y) return true;
            else if (v1.uv.y > v2.uv.y) return false;
            return false;
        }
    };

    void createMesh(
            // in parameters
            const std::string &name,
            const std::vector<Point> &vertices,
            const std::vector<Normal> &normals,
            const std::vector<Point2> &texcoords,
            const std::vector<OBJTriangle> &triangles,
            const std::string &materialName,
            const Transform &objectToWorld,
            std::vector<Vertex> &vertexBuffer,
            // out parameters
            std::vector<TriMesh *> &m_meshes,
            std::vector<std::string> &m_materialAssignment) {
        if (triangles.size() == 0)
            return;
        typedef std::map<Vertex, uint32_t, vertex_key_order> VertexMapType;
        VertexMapType vertexMap;

        vertexBuffer.reserve(vertices.size());
        size_t numMerged = 0;
        AABB aabb;
        bool hasTexcoords = false;
        bool hasNormals = false;
        ref<Timer> timer = new Timer();
        vertexBuffer.clear();

        /* Collapse the mesh into a more usable form */
        Triangle *triangleArray = new Triangle[triangles.size()];
        for (uint32_t i=0; i<triangles.size(); i++) {
            Triangle tri;
            for (uint32_t j=0; j<3; j++) {
                int vertexId = triangles[i].p[j];
                int normalId = triangles[i].n[j];
                int uvId = triangles[i].uv[j];
                uint32_t key;

                Vertex vertex;
                if (vertexId < 0)
                    vertexId += (int) vertices.size() + 1;
                if (normalId < 0)
                    normalId += (int) normals.size() + 1;
                if (uvId < 0)
                    uvId += (int) texcoords.size() + 1;

                if (vertexId > (int) vertices.size() || vertexId <= 0)
                    Log(EError, "Out of bounds: tried to access vertex %i (max: %i)", vertexId, (int) vertices.size());

                vertex.p = objectToWorld(vertices[vertexId-1]);
                aabb.expandBy(vertex.p);

                if (normalId != 0) {
                    if (normalId > (int) normals.size() || normalId < 0)
                        Log(EError, "Out of bounds: tried to access normal %i (max: %i)", normalId, (int) normals.size());
                    vertex.n = objectToWorld(normals[normalId-1]);
                    if (!vertex.n.isZero())
                        vertex.n = normalize(vertex.n);
                    hasNormals = true;
                } else {
                    vertex.n = Normal(0.0f);
                }

                if (uvId != 0) {
                    if (uvId > (int) texcoords.size() || uvId < 0)
                        Log(EError, "Out of bounds: tried to access uv %i (max: %i)", uvId, (int) texcoords.size());
                    vertex.uv = texcoords[uvId-1];
                    hasTexcoords = true;
                } else {
                    vertex.uv = Point2(0.0f);
                }

                VertexMapType::iterator it = vertexMap.find(vertex);
                if (it != vertexMap.end()) {
                    key = it->second;
                    numMerged++;
                } else {
                    key = (uint32_t) vertexBuffer.size();
                    vertexMap[vertex] = key;
                    vertexBuffer.push_back(vertex);
                }

                tri.idx[j] = key;
            }
            triangleArray[i] = tri;
        }

        ref<TriMesh> mesh = new TriMesh(name,
            triangles.size(), vertexBuffer.size(),
            hasNormals, hasTexcoords, false,
            m_flipNormals, m_faceNormals);

        std::copy(triangleArray, triangleArray+triangles.size(), mesh->getTriangles());

        Point    *target_positions = mesh->getVertexPositions();
        Normal   *target_normals   = mesh->getVertexNormals();
        Point2   *target_texcoords = mesh->getVertexTexcoords();

        mesh->getAABB() = aabb;

        for (size_t i=0; i<vertexBuffer.size(); i++) {
            *target_positions++ = vertexBuffer[i].p;
            if (hasNormals)
                *target_normals++ = vertexBuffer[i].n;
            if (hasTexcoords)
                *target_texcoords++ = vertexBuffer[i].uv;
        }

        mesh->incRef();
        m_materialAssignment.push_back(materialName);
        m_meshes.push_back(mesh);
        Log(EInfo, "%s: " SIZE_T_FMT " triangles, " SIZE_T_FMT
            " vertices (merged " SIZE_T_FMT " vertices).", name.c_str(),
            triangles.size(), vertexBuffer.size(), numMerged);
    }

    void loadWavefrontOBJ(
            // out parameters
            std::vector<TriMesh *> &m_meshes,
            std::vector<std::string> &m_materialAssignment,
            // in parameters
            std::string m_name,
            fs::ifstream &is,
            ref<FileResolver> &fileResolver,
            Transform &objectToWorld,
            bool m_collapse = false,
            bool flipTexCoords = true,
            bool loadMaterials = true,
            int shapeIndex = -1)
    {

        std::string buf;
        std::vector<Point> vertices;
        std::vector<Normal> normals;
        std::vector<Point2> texcoords;
        std::vector<OBJTriangle> triangles;
        std::string name = m_name, line;
        std::set<std::string> geomNames;
        std::vector<Vertex> vertexBuffer;
        fs::path materialLibrary;
        int geomIndex = 0;
        bool nameBeforeGeometry = false;
        std::string materialName;

        while (is.good() && !is.eof() && fetch_line(is, line)) {
            std::istringstream iss(line);
            if (!(iss >> buf))
                continue;

            if (buf == "v") {
                /* Parse + transform vertices */
                Point p;
                iss >> p.x >> p.y >> p.z;
                vertices.push_back(p);
            } else if (buf == "vn") {
                Normal n;
                iss >> n.x >> n.y >> n.z;
                normals.push_back(n);
            } else if (buf == "g" && !m_collapse) {
                std::string targetName;
                std::string newName = trim(line.substr(1, line.length()-1));

                /* There appear to be two different conventions
                   for specifying object names in OBJ file -- try
                   to detect which one is being used */
                if (nameBeforeGeometry)
                    // Save geometry under the previously specified name
                    targetName = name;
                else
                    targetName = newName;

                if (triangles.size() > 0) {
                    /// make sure that we have unique names
                    if (geomNames.find(targetName) != geomNames.end())
                        targetName = formatString("%s_%i", targetName.c_str(), geomIndex);
                    geomIndex += 1;
                    geomNames.insert(targetName);
                    if (shapeIndex < 0 || geomIndex-1 == shapeIndex)
                        createMesh(targetName, vertices, normals, texcoords,
                            triangles, materialName, objectToWorld, vertexBuffer,
                            m_meshes, m_materialAssignment);
                    triangles.clear();
                } else {
                    nameBeforeGeometry = true;
                }
                name = newName;
            } else if (buf == "usemtl") {
                /* Flush if necessary */
                if (triangles.size() > 0 && !m_collapse) {
                    /// make sure that we have unique names
                    if (geomNames.find(name) != geomNames.end())
                        name = formatString("%s_%i", name.c_str(), geomIndex);
                    geomIndex += 1;
                    geomNames.insert(name);
                    if (shapeIndex < 0 || geomIndex-1 == shapeIndex)
                        createMesh(name, vertices, normals, texcoords,
                            triangles, materialName, objectToWorld, vertexBuffer,
                            m_meshes, m_materialAssignment);
                    triangles.clear();
                    name = m_name;
                }

                materialName = trim(line.substr(6, line.length()-1));
            } else if (buf == "mtllib") {
                materialLibrary = fileResolver->resolve(trim(line.substr(6, line.length()-1)));
            } else if (buf == "vt") {
                Float u, v;
                iss >> u >> v;
                if (flipTexCoords)
                    v = 1-v;
                texcoords.push_back(Point2(u, v));
            } else if (buf == "f") {
                std::string  tmp;
                OBJTriangle t;
                iss >> tmp; parse(t, 0, tmp);
                iss >> tmp; parse(t, 1, tmp);
                iss >> tmp; parse(t, 2, tmp);
                triangles.push_back(t);
                /* Handle n-gons assuming a convex shape */
                while (iss >> tmp) {
                    t.p[1] = t.p[2];
                    t.uv[1] = t.uv[2];
                    t.n[1] = t.n[2];
                    parse(t, 2, tmp);
                    triangles.push_back(t);
                }
            } else {
                /* Ignore */
            }
        }
        if (geomNames.find(name) != geomNames.end())
            /// make sure that we have unique names
            name = formatString("%s_%i", m_name.c_str(), geomIndex);

        if (shapeIndex < 0 || geomIndex-1 == shapeIndex)
            createMesh(name, vertices, normals, texcoords,
                triangles, materialName, objectToWorld, vertexBuffer,
                m_meshes, m_materialAssignment);

        if (!materialLibrary.empty() && loadMaterials)
            loadMaterialLibrary(fileResolver, materialLibrary);
    }


    void serialize(Stream *stream, InstanceManager *manager) const {
        Shape::serialize(stream, manager);

        m_aabb.serialize(stream);
        stream->writeUInt((uint32_t) m_meshes.size());
        for (size_t i=0; i<m_meshes.size(); ++i)
            manager->serialize(stream, m_meshes[i]);
    }

    virtual ~GrainMesh() {
        for (size_t i=0; i<m_meshes.size(); ++i)
            m_meshes[i]->decRef();   
    }


    void addChild(const std::string &name, ConfigurableObject *child) {
        addChild(name, child, true);
    }

    void addChild(const std::string &name, ConfigurableObject *child, bool warn) {
        const Class *cClass = child->getClass();
        if (cClass->derivesFrom(MTS_CLASS(BSDF))) {
            Shape::addChild(name, child);
            Assert(m_meshes.size() > 0);
            if (name == "") {
                for (size_t i=0; i<m_meshes.size(); ++i)
                    m_meshes[i]->addChild(name, child);
            } else {
                bool found = false;
                for (size_t i=0; i<m_meshes.size(); ++i) {
                    if (m_materialAssignment[i] == name) {
                        found = true;
                        m_meshes[i]->addChild(name, child);
                    }
                }
                if (!found && warn)
                    Log(EWarn, "Attempted to register the material named "
                        "'%s', which does not occur in the OBJ file!", name.c_str());
            }
            m_bsdf->setParent(NULL);
        } else if (cClass->derivesFrom(MTS_CLASS(Emitter))) {
            if (m_meshes.size() > 1)
                ; //TODO any side-effects?
                //Log(EError, "Cannot attach an emitter to an OBJ file "
                    //"containing multiple objects!");
            m_emitter = static_cast<Emitter *>(child);
            child->setParent(m_meshes[0]);
            m_meshes[0]->addChild(name, child);
        } else if (cClass->derivesFrom(MTS_CLASS(Sensor))) {
            if (m_meshes.size() > 1)
                ; //TODO any side-effects?
                //Log(EError, "Cannot attach an sensor to an OBJ file "
                    //"containing multiple objects!");
            m_sensor = static_cast<Sensor *>(child);
            child->setParent(m_meshes[0]);
            m_meshes[0]->addChild(name, child);
        } else if (cClass->derivesFrom(MTS_CLASS(Subsurface))) {
            Assert(m_subsurface == NULL);
            m_subsurface = static_cast<Subsurface *>(child);
            for (size_t i=0; i<m_meshes.size(); ++i) {
                child->setParent(m_meshes[i]);
                m_meshes[i]->addChild(name, child);
            }
        } else if (cClass->derivesFrom(MTS_CLASS(Medium))) {
            Shape::addChild(name, child);
            for (size_t i=0; i<m_meshes.size(); ++i)
                m_meshes[i]->addChild(name, child);
        } else {
            Shape::addChild(name, child);
        }
    }

#if 0
    bool isCompound() const {
        return true;
    }
#endif

    Shape *getElement(int index) {
        if (index >= (int) m_meshes.size())
            return NULL;
        Shape *shape = m_meshes[index];
        BSDF *bsdf = shape->getBSDF();
        Emitter *emitter = shape->getEmitter();
        Subsurface *subsurface = shape->getSubsurface();
        if (bsdf)
            bsdf->setParent(shape);
        if (emitter)
            emitter->setParent(shape);
        if (subsurface)
            subsurface->setParent(shape);
        return shape;
    }

    AABB getAABB() const {
        return m_volume_aabb;
        //return m_aabb;
    }

    Float getSurfaceArea() const {
        Float sa = 0;
        for (size_t i=0; i<m_meshes.size(); ++i)
            sa += m_meshes[i]->getSurfaceArea();
        return sa;
    }

    size_t getPrimitiveCount() const {
        return 1; //ying

        size_t result = 0;
        for (size_t i=0; i<m_meshes.size(); ++i)
            result += m_meshes[i]->getPrimitiveCount();
        return result;
    }

    size_t getEffectivePrimitiveCount() const {
        return 1; //ying

        size_t result = 0;
        for (size_t i=0; i<m_meshes.size(); ++i)
            result += m_meshes[i]->getEffectivePrimitiveCount();
        return result;
    }


    MTS_DECLARE_CLASS()
private:
    std::vector<TriMesh *> m_meshes;
    std::vector<std::string> m_materialAssignment;
    bool m_flipNormals, m_faceNormals;
    AABB m_aabb, m_volume_aabb;
    bool m_collapse;
};
#endif


# if TEST
class GrainMesh : public Shape {
private:
    Transform m_objectToWorld;
    Transform m_worldToObject;
    Float m_radius, m_length, m_invSurfaceArea;
    bool m_flipNormals;

public:
    GrainMesh(const Properties &props) : Shape(props) {
# if 0
        // int dimension = props.getInteger("dimension", 128);

        ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver()->clone();
        fs::path path = fileResolver->resolve(props.getString("filename"));

        m_name = path.stem().string();

        /* Object-space -> World-space transformation */
        //Transform objectToWorld = props.getTransform("toWorld", Transform());

        /* Load the geometry */
        Log(EInfo, "Loading binvox from \"%s\" ..", path.filename().string().c_str());

        fileResolver->prependPath(fs::absolute(path).parent_path());

        ref<Timer> timer = new Timer();

        read_binvox(path.string()+".binvox");
        read_binvox(path.string()+".hollow.binvox", false);

        Log(EInfo, "Done with \"%s\" (took %i ms)", path.filename().string().c_str(), timer->getMilliseconds());

        m_obj = new WavefrontOBJ(props);
        AABB aabb = m_obj->getAABB();
        Log(EInfo, "Mesh aabb: %f %f %f, %f %f %f",
            aabb.min.x, aabb.min.y, aabb.min.z, aabb.max.x, aabb.max.y, aabb.max.z);
#endif

        Float radius = props.getFloat("radius", 1.0f);
        Point p1 = props.getPoint("p0", Point(0.0f, 0.0f, 0.0f));
        Point p2 = props.getPoint("p1", Point(0.0f, 0.0f, 1.0f));
        Vector d = p2 - p1;
        Float length = d.length();
        m_objectToWorld =
            Transform::translate(Vector(p1)) *
            Transform::fromFrame(Frame(d / length)) *
            Transform::scale(Vector(radius, radius, length));

        if (props.hasProperty("toWorld"))
            m_objectToWorld = props.getTransform("toWorld") * m_objectToWorld;

        /// Are the cylinder normals pointing inwards? default: no
        m_flipNormals = props.getBoolean("flipNormals", false);

        // Remove the scale from the object-to-world transform
        m_radius = m_objectToWorld(Vector(1,0,0)).length();
        m_length = m_objectToWorld(Vector(0,0,1)).length();
        m_objectToWorld = m_objectToWorld * Transform::scale(
            Vector(1/m_radius, 1/m_radius, 1/m_length));

        m_worldToObject = m_objectToWorld.inverse();
        m_invSurfaceArea = 1/(2*M_PI*m_radius*m_length);
        Assert(m_length > 0 && m_radius > 0);
    }

    GrainMesh(Stream *stream, InstanceManager *manager)
        : Shape(stream, manager) {
        m_objectToWorld = Transform(stream);
        m_radius = stream->readFloat();
        m_length = stream->readFloat();
        m_flipNormals = stream->readBool();
        m_worldToObject = m_objectToWorld.inverse();
        m_invSurfaceArea = 1/(2*M_PI*m_radius*m_length);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        Shape::serialize(stream, manager);
        m_objectToWorld.serialize(stream);
        stream->writeFloat(m_radius);
        stream->writeFloat(m_length);
        stream->writeBool(m_flipNormals);
    }

    bool rayIntersect(const Ray &_ray, Float mint, Float maxt, Float &t, void *temp) const {
        Ray ray;
        /* Transform into the local coordinate system and normalize */
        m_worldToObject(_ray, ray);

        const double
            ox = ray.o.x,
            oy = ray.o.y,
            dx = ray.d.x,
            dy = ray.d.y;

        const double A = dx*dx + dy*dy;
        const double B = 2 * (dx*ox + dy*oy);
        const double C = ox*ox + oy*oy - m_radius*m_radius;

        double nearT, farT;
        if (!solveQuadraticDouble(A, B, C, nearT, farT))
            return false;

        if (!(nearT <= maxt && farT >= mint)) /* NaN-aware conditionals */
            return false;

        const double zPosNear = ray.o.z + ray.d.z * nearT;
        const double zPosFar = ray.o.z + ray.d.z * farT;

        if (zPosNear >= 0 && zPosNear <= m_length && nearT >= mint) {
            t = (Float) nearT;
        } else if (zPosFar >= 0 && zPosFar <= m_length) {
            if (farT > maxt)
                return false;
            t = (Float) farT;
        } else {
            return false;
        }

        return true;
    }

    bool rayIntersect(const Ray &_ray, Float mint, Float maxt) const {
        Ray ray;
        /* Transform into the local coordinate system and normalize */
        m_worldToObject(_ray, ray);

        const double
            ox = ray.o.x,
            oy = ray.o.y,
            dx = ray.d.x,
            dy = ray.d.y;

        const double A = dx*dx + dy*dy;
        const double B = 2 * (dx*ox + dy*oy);
        const double C = ox*ox + oy*oy - m_radius*m_radius;

        double nearT, farT;
        if (!solveQuadraticDouble(A, B, C, nearT, farT))
            return false;

        if (nearT > maxt || farT < mint)
            return false;

        const double zPosNear = ray.o.z + ray.d.z * nearT;
        const double zPosFar = ray.o.z + ray.d.z * farT;
        if (zPosNear >= 0 && zPosNear <= m_length && nearT >= mint) {
            return true;
        } else if (zPosFar >= 0 && zPosFar <= m_length && farT <= maxt) {
            return true;
        } else {
            return false;
        }
    }

    void fillIntersectionRecord(const Ray &ray,
            const void *temp, Intersection &its) const {
        its.p = ray(its.t);
        Point local = m_worldToObject(its.p);

        Float phi = std::atan2(local.y, local.x);
        if (phi < 0)
            phi += 2*M_PI;
        its.uv.x = phi / (2*M_PI);
        its.uv.y = local.z / m_length;

        Vector dpdu = Vector(-local.y, local.x, 0) * (2*M_PI);
        Vector dpdv = Vector(0, 0, m_length);
        its.shape = this;
        its.dpdu = m_objectToWorld(dpdu);
        its.dpdv = m_objectToWorld(dpdv);
        its.geoFrame.s = normalize(its.dpdu);
        its.geoFrame.t = normalize(its.dpdv);
        its.geoFrame.n = Normal(cross(its.geoFrame.s, its.geoFrame.t));

        /* Migitate roundoff error issues by a normal shift of the computed intersection point */
        its.p += its.geoFrame.n * (m_radius - std::sqrt(local.x*local.x+local.y*local.y));

        if (m_flipNormals)
            its.geoFrame.n *= -1;
        its.shFrame.n = its.geoFrame.n;
        its.hasUVPartials = false;
        its.instance = NULL;
        its.time = ray.time;
    }

    void samplePosition(PositionSamplingRecord &pRec, const Point2 &sample) const {
        Float sinTheta, cosTheta;
        math::sincos(sample.y * (2 * M_PI), &sinTheta, &cosTheta);

        Point p(cosTheta*m_radius, sinTheta*m_radius, sample.x * m_length);
        Normal n(cosTheta, sinTheta, 0.0f);

        if (m_flipNormals)
            n *= -1;

        pRec.p = m_objectToWorld(p);
        pRec.n = normalize(m_objectToWorld(n));
        pRec.pdf = m_invSurfaceArea;
        pRec.measure = EArea;
    }

    Float pdfPosition(const PositionSamplingRecord &pRec) const {
        return m_invSurfaceArea;
    }

    inline AABB getAABB() const {
        Vector x1 = m_objectToWorld(Vector(m_radius, 0, 0));
        Vector x2 = m_objectToWorld(Vector(0, m_radius, 0));
        Point p0 = m_objectToWorld(Point(0, 0, 0));
        Point p1 = m_objectToWorld(Point(0, 0, m_length));
        AABB result;

        /* To bound the cylinder, it is sufficient to find the
           smallest box containing the two circles at the endpoints.
           This can be done component-wise as follows */

        for (int i=0; i<3; ++i) {
            Float range = std::sqrt(x1[i]*x1[i] + x2[i]*x2[i]);

            result.min[i] = std::min(std::min(result.min[i],
                        p0[i]-range), p1[i]-range);
            result.max[i] = std::max(std::max(result.max[i],
                        p0[i]+range), p1[i]+range);
        }

        return result;
    }

    /**
     * Compute the ellipse created by the intersection of an infinite
     * cylinder and a plane. Returns false in the degenerate case.
     * Based on:
     * www.geometrictools.com/Documentation/IntersectionCylinderPlane.pdf
     */
    bool intersectCylPlane(Point planePt, Normal planeNrml,
            Point cylPt, Vector cylD, Float radius, Point &center,
            Vector *axes, Float *lengths) const {
        if (absDot(planeNrml, cylD) < Epsilon)
            return false;

        Vector B, A = cylD - dot(cylD, planeNrml)*planeNrml;

        Float length = A.length();
        if (length != 0) {
            A /= length;
            B = cross(planeNrml, A);
        } else {
            coordinateSystem(planeNrml, A, B);
        }

        Vector delta = planePt - cylPt,
               deltaProj = delta - cylD*dot(delta, cylD);

        Float aDotD = dot(A, cylD);
        Float bDotD = dot(B, cylD);
        Float c0 = 1-aDotD*aDotD;
        Float c1 = 1-bDotD*bDotD;
        Float c2 = 2*dot(A, deltaProj);
        Float c3 = 2*dot(B, deltaProj);
        Float c4 = dot(delta, deltaProj) - radius*radius;

        Float lambda = (c2*c2/(4*c0) + c3*c3/(4*c1) - c4)/(c0*c1);

        Float alpha0 = -c2/(2*c0),
              beta0 = -c3/(2*c1);

        lengths[0] = std::sqrt(c1*lambda),
        lengths[1] = std::sqrt(c0*lambda);

        center = planePt + alpha0 * A + beta0 * B;
        axes[0] = A;
        axes[1] = B;
        return true;
    }

    AABB intersectCylFace(int axis,
            const Point &min, const Point &max,
            const Point &cylPt, const Vector &cylD) const {
        int axis1 = (axis + 1) % 3;
        int axis2 = (axis + 2) % 3;

        Normal planeNrml(0.0f);
        planeNrml[axis] = 1;

        Point ellipseCenter;
        Vector ellipseAxes[2];
        Float ellipseLengths[2];

        AABB aabb;
        if (!intersectCylPlane(min, planeNrml, cylPt, cylD, m_radius,
            ellipseCenter, ellipseAxes, ellipseLengths)) {
            /* Degenerate case -- return an invalid AABB. This is
               not a problem, since one of the other faces will provide
               enough information to arrive at a correct clipped AABB */
            return aabb;
        }

        /* Intersect the ellipse against the sides of the AABB face */
        for (int i=0; i<4; ++i) {
            Point p1, p2;
            p1[axis] = p2[axis] = min[axis];
            p1[axis1] = ((i+1) & 2) ? min[axis1] : max[axis1];
            p1[axis2] = ((i+0) & 2) ? min[axis2] : max[axis2];
            p2[axis1] = ((i+2) & 2) ? min[axis1] : max[axis1];
            p2[axis2] = ((i+1) & 2) ? min[axis2] : max[axis2];

            Point2 p1l(
                dot(p1 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
                dot(p1 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);
            Point2 p2l(
                dot(p2 - ellipseCenter, ellipseAxes[0]) / ellipseLengths[0],
                dot(p2 - ellipseCenter, ellipseAxes[1]) / ellipseLengths[1]);

            Vector2 rel = p2l-p1l;
            Float A = dot(rel, rel);
            Float B = 2*dot(Vector2(p1l), rel);
            Float C = dot(Vector2(p1l), Vector2(p1l))-1;

            Float x0, x1;
            if (solveQuadratic(A, B, C, x0, x1)) {
                if (x0 >= 0 && x0 <= 1)
                    aabb.expandBy(p1+(p2-p1)*x0);
                if (x1 >= 0 && x1 <= 1)
                    aabb.expandBy(p1+(p2-p1)*x1);
            }
        }

        ellipseAxes[0] *= ellipseLengths[0];
        ellipseAxes[1] *= ellipseLengths[1];
        AABB faceBounds(min, max);

        /* Find the componentwise maxima of the ellipse */
        for (int i=0; i<2; ++i) {
            int j = (i==0) ? axis1 : axis2;
            Float alpha = ellipseAxes[0][j], beta = ellipseAxes[1][j];
            Float tmp = 1 / std::sqrt(alpha*alpha + beta*beta);
            Float cosTheta = alpha * tmp, sinTheta = beta*tmp;

            Point p1 = ellipseCenter + cosTheta*ellipseAxes[0] + sinTheta*ellipseAxes[1];
            Point p2 = ellipseCenter - cosTheta*ellipseAxes[0] - sinTheta*ellipseAxes[1];

            if (faceBounds.contains(p1))
                aabb.expandBy(p1);
            if (faceBounds.contains(p2))
                aabb.expandBy(p2);
        }

        return aabb;
    }



    AABB getClippedAABB(const AABB &box) const {
        /* Compute a base bounding box */
        AABB base(getAABB());
        base.clip(box);

        Point cylPt = m_objectToWorld(Point(0, 0, 0));
        Vector cylD(m_objectToWorld(Vector(0, 0, 1)));

        /* Now forget about the cylinder ends and
           intersect an infinite cylinder with each AABB face */
        AABB clippedAABB;
        clippedAABB.expandBy(intersectCylFace(0,
                Point(base.min.x, base.min.y, base.min.z),
                Point(base.min.x, base.max.y, base.max.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(0,
                Point(base.max.x, base.min.y, base.min.z),
                Point(base.max.x, base.max.y, base.max.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(1,
                Point(base.min.x, base.min.y, base.min.z),
                Point(base.max.x, base.min.y, base.max.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(1,
                Point(base.min.x, base.max.y, base.min.z),
                Point(base.max.x, base.max.y, base.max.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(2,
                Point(base.min.x, base.min.y, base.min.z),
                Point(base.max.x, base.max.y, base.min.z),
                cylPt, cylD));

        clippedAABB.expandBy(intersectCylFace(2,
                Point(base.min.x, base.min.y, base.max.z),
                Point(base.max.x, base.max.y, base.max.z),
                cylPt, cylD));

        clippedAABB.clip(box);
        return clippedAABB;
    }

    ref<TriMesh> createTriMesh() {
        /// Choice of discretization
        const size_t phiSteps = 20;
        const Float dPhi   = (2*M_PI) / phiSteps;

        ref<TriMesh> mesh = new TriMesh("Cylinder approximation",
            phiSteps*2, phiSteps*2, true, false, false);

        Point *vertices = mesh->getVertexPositions();
        Normal *normals = mesh->getVertexNormals();
        Triangle *triangles = mesh->getTriangles();
        size_t triangleIdx = 0, vertexIdx = 0;

        for (size_t phi=0; phi<phiSteps; ++phi) {
            Float sinPhi = std::sin(phi * dPhi);
            Float cosPhi = std::cos(phi * dPhi);
            uint32_t idx0 = (uint32_t) vertexIdx, idx1 = idx0+1;
            uint32_t idx2 = (vertexIdx+2) % (2*phiSteps), idx3 = idx2+1;
            normals[vertexIdx] = m_objectToWorld(Normal(cosPhi, sinPhi, 0) * (m_flipNormals ? (Float) -1 : (Float) 1));
            vertices[vertexIdx++] = m_objectToWorld(Point(cosPhi*m_radius, sinPhi*m_radius, 0));
            normals[vertexIdx] = m_objectToWorld(Normal(cosPhi, sinPhi, 0) * (m_flipNormals ? (Float) -1 : (Float) 1));
            vertices[vertexIdx++] = m_objectToWorld(Point(cosPhi*m_radius, sinPhi*m_radius, m_length));

            triangles[triangleIdx].idx[0] = idx0;
            triangles[triangleIdx].idx[1] = idx2;
            triangles[triangleIdx].idx[2] = idx1;
            triangleIdx++;
            triangles[triangleIdx].idx[0] = idx1;
            triangles[triangleIdx].idx[1] = idx2;
            triangles[triangleIdx].idx[2] = idx3;
            triangleIdx++;
        }

        mesh->copyAttachments(this);
        mesh->configure();

        return mesh.get();
    }

#if 0
    AABB getAABB() const {
        const Point a = m_objectToWorld(Point(0, 0, 0));
        const Point b = m_objectToWorld(Point(0, 0, m_length));

        const Float r = m_radius;
        AABB result;
        result.expandBy(a - Vector(r, r, r));
        result.expandBy(a + Vector(r, r, r));
        result.expandBy(b - Vector(r, r, r));
        result.expandBy(b + Vector(r, r, r));
        return result;
    }
#endif

    Float getSurfaceArea() const {
        return 2*M_PI*m_radius*m_length;
    }

    void getNormalDerivative(const Intersection &its,
            Vector &dndu, Vector &dndv, bool shadingFrame) const {
        dndu = its.dpdu / (m_radius * (m_flipNormals ? -1 : 1));
        dndv = Vector(0.0f);
    }

    size_t getPrimitiveCount() const {
        return 1;
    }

    size_t getEffectivePrimitiveCount() const {
        return 1;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "GrainMesh[" << endl
            << "  radius = " << m_radius << "," << endl
            << "  length = " << m_length << "," << endl
            << "  objectToWorld = " << indent(m_objectToWorld.toString()) << "," << endl
            << "  bsdf = " << indent(m_bsdf.toString()) << "," << endl;
        if (isMediumTransition())
            oss << "  interiorMedium = " << indent(m_interiorMedium.toString()) << "," << endl
                << "  exteriorMedium = " << indent(m_exteriorMedium.toString()) << "," << endl;
        oss << "  emitter = " << indent(m_emitter.toString()) << "," << endl
            << "  sensor = " << indent(m_sensor.toString()) << "," << endl
            << "  subsurface = " << indent(m_subsurface.toString())
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};
#endif

MTS_IMPLEMENT_CLASS_S(GrainMesh, false, Shape)
MTS_EXPORT_PLUGIN(GrainMesh, "GrainMesh intersection primitive");
MTS_NAMESPACE_END

