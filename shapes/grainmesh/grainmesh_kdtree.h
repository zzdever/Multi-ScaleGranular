#ifndef __GRAINMESH_KDTREE__
#define __GRAINMESH_KDTREE__

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
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/aabb.h>

#include "tempinfo.h"

using namespace mitsuba;

namespace grain {

struct KDTriangle {
    //uint32_t idx[3];
    Point m_points[3];
    bool m_hasNormal;
    Normal m_normals[3];
    bool m_hasUV;
    Point2 m_uvs[3];

    KDTriangle(Point p1, Point p2, Point p3) {
        m_hasNormal = false;
        m_hasUV = false;
        m_points[0] = p1;
        m_points[1] = p2;
        m_points[2] = p3;
    }

    inline void addNormal(Normal n1, Normal n2, Normal n3) {
        m_hasUV = true;
        m_normals[0] = n1;
        m_normals[1] = n2;
        m_normals[2] = n3;
    }

    inline void addUV(Point2 uv1, Point2 uv2, Point2 uv3) {
        m_hasUV = true;
        m_uvs[0] = uv1;
        m_uvs[1] = uv2;
        m_uvs[2] = uv3;
    }

    inline AABB getAABB() const {
        AABB result(m_points[0]);
        result.expandBy(m_points[1]);
        result.expandBy(m_points[2]);
        return result;
    }

    inline Point get_midpoint() const {
        return (m_points[0] + m_points[1] + m_points[2]) / 3.0;
    }

    inline Normal getNormal() const {
        return (m_normals[0] + m_normals[1] + m_normals[2]) / 3.0;
    }

    inline Point2 getUV() const {
        return (m_uvs[0] + m_uvs[1] + m_uvs[2]) / 3.0;
    }

    inline TangentSpace getUVTangent() const {
        // TODO use a correct one
        TangentSpace t;
        t.dpdu = (m_points[1] - m_points[0]);
        t.dpdv = (m_points[2] - m_points[0]);
        return t;
    }

    inline Vector operator=(Point p) {
        return Vector(p.x, p.y, p.z);
    }

    bool hit(const Ray& ray, float& to, float& tmin) {
        // return 1==hit_(ray, to, tmin);
        Float u,v,t;
        bool res = rayIntersect(m_points[0], m_points[1], m_points[2],
                ray, u, v, t) ;
        if(res && t>=tmin) {
            to = (float)t;
            return true;
        }
        return false;
    }


#if 0
#define SMALL_NUM   0.00000001 // anything that avoids division overflow

    //    Return: -1 = triangle is degenerate (a segment or point)
    //             0 =  disjoint (no intersect)
    //             1 =  intersect in unique point
    //             2 =  are in the same plane
    inline int hit_(const Ray& ray, float& ot, float& tmin)
    {
        Vector    u, v, n;              // triangle vectors
        Vector    dir, w0, w;           // ray vectors
        float     r, a, b;              // params to calc ray-plane intersect

        // get triangle edge vectors and plane normal
        //u = T.V1 - T.V0;
        //v = T.V2 - T.V0;
        u = m_points[1] - m_points[0];
        v = m_points[2] - m_points[0];
        n = cross(u, v);              // cross product
        if (n == Vector(0,0,0))             // triangle is degenerate
            return -1;                  // do not deal with this case

        dir = ray.d;              // ray direction vector
        w0 = ray.o - m_points[0];
        a = -dot(n, w0);
        b = dot(n, dir);
        if (fabs(b) < SMALL_NUM) {     // ray is  parallel to triangle plane
            if (a == 0)                 // ray lies in triangle plane
                return 2;
            else return 0;              // ray disjoint from plane
        }

        // get intersect point of ray with triangle plane
        r = a / b;
        if (r < 0.0)                    // ray goes away from triangle
            return 0;                   // => no intersect
        // for a segment, also test if (r > 1.0) => no intersect

        //*I = R.P0 + r * dir;            // intersect point of ray and plane
        ot = r;
        if(ot <= tmin)
            return 0; // the tri is behind the ray

        // is I inside T?
        float uu, uv, vv, wu, wv, D;
        uu = dot(u, u);
        uv = dot(u, v);
        vv = dot(v, v);
        //w = *I - T.V0;
        w = ray.o + r * dir - m_points[0];
        wu = dot(w, u);
        wv = dot(w, v);
        D = uv * uv - uu * vv;

        // get and test parametric coords
        float s, t;
        s = (uv * wv - vv * wu) / D;
        if (s < 0.0 || s > 1.0)         // I is outside T
            return 0;
        t = (uv * wu - uu * wv) / D;
        if (t < 0.0 || (s + t) > 1.0)  // I is outside T
            return 0;

        return 1;                       // I is in T
    }
#endif

#if 0
    /**
     * \brief Returns the axis-aligned bounding box of a triangle after it has
     * clipped to the extends of another given AABB.
     *
     * This function uses the Sutherland-Hodgman algorithm to calculate the
     * convex polygon that is created when applying all 6 AABB splitting
     * planes to the triangle. Afterwards, the AABB of the newly created
     * convex polygon is returned. This function is an important component
     * for efficiently creating 'Perfect Split' KD-trees. For more detail,
     * see "On building fast kd-Trees for Ray Tracing, and on doing
     * that in O(N log N)" by Ingo Wald and Vlastimil Havran
     */
    AABB getClippedAABB(const Point *positions, const AABB &aabb) const;

    // Returns the bounding sphere of the triangle
    inline BSphere getBSphere(const Point *positions) const {
        Vector a = (positions[idx[1]] - positions[idx[0]]);
        Vector b = (positions[idx[2]] - positions[idx[0]]);
        Float a2 = dot(a, a);
        Float b2 = dot(b, b);
        Float da = std::sqrt(a2);
        Float db = std::sqrt(b2);
        Vector axb = cross(a, b);
        Float axb2 = dot(axb, axb);
        Float daxb = std::sqrt(axb2);
        return BSphere(positions[idx[0]] + cross(a2 * b - b2 * a, axb) / (2 * axb2),
                       da * db * (a - b).length() / (2 * daxb));
    }

    /// Uniformly sample a point on the triangle and return its normal and UV coordinates
    Point sample(const Point *positions, const Normal *normals,
            const Point2 *texCoords, Normal &n, Point2 &uv,
            const Point2 &seed) const;

    /// Calculate the surface area of this triangle
    Float surfaceArea(const Point *positions) const;

#endif
    /** \brief Ray-triangle intersection test
     *
     * Uses the algorithm by Moeller and Trumbore discussed at
     * <tt>http://www.acm.org/jgt/papers/MollerTrumbore97/code.html</tt>.
     *
     * \param p0
     *    Position of the first vertex
     * \param p1
     *    Position of the second vertex
     * \param p2
     *    Position of the third vertex
     * \param ray
     *    The ray segment to be used for the intersection query
     * \param t
     *    Upon success, \a t contains the distance from the ray origin to the
     *    intersection point,
     * \param u
     *   Upon success, \c u will contain the 'U' component of the intersection
     *   in barycentric coordinates
     * \param v
     *   Upon success, \c v will contain the 'V' component of the intersection
     *   in barycentric coordinates
     * \return
     *   \c true if an intersection has been detected
     */
    FINLINE static bool rayIntersect(const Point &p0, const Point &p1, const Point &p2,
        const Ray &ray, Float &u, Float &v, Float &t) {
        /* Find vectors for two edges sharing */
        Vector edge1 = p1 - p0, edge2 = p2 - p0;

        /* Begin calculating determinant - also used to calculate U parameter */
        Vector pvec = cross(ray.d, edge2);

        Float det = dot(edge1, pvec);
        if (det == 0)
            return false;
        Float inv_det = 1.0f / det;

        /* Calculate distance from v[0] to ray origin */
        Vector tvec = ray.o - p0;

        /* Calculate U parameter and test bounds */
        u = dot(tvec, pvec) * inv_det;
        if (u < 0.0 || u > 1.0)
            return false;

        /* Prepare to test V parameter */
        Vector qvec = cross(tvec, edge1);

        /* Calculate V parameter and test bounds */
        v = dot(ray.d, qvec) * inv_det;

        /* Inverted comparison (to catch NaNs) */
        if (v >= 0.0 && u + v <= 1.0) {
            /* ray intersects triangle -> compute t */
            t = dot(edge2, qvec) * inv_det;

            return true;
        }

        return false;
    }

#if 0
    /** \brief Ray-triangle intersection test
     *
     * Uses the algorithm by Moeller and Trumbore discussed at
     * <tt>http://www.acm.org/jgt/papers/MollerTrumbore97/code.html</tt>.
     *
     * \param positions
     *    Pointer to the vertex positions of the underlying triangle mesh
     * \param index
     *    Index of the triangle that should be intersected
     * \param ray
     *    The ray segment to be used for the intersection query
     * \param t
     *    Upon success, \a t contains the distance from the ray origin to the
     *    intersection point,
     * \param u
     *   Upon success, \c u will contain the 'U' component of the intersection
     *   in barycentric coordinates
     * \param v
     *   Upon success, \c v will contain the 'V' component of the intersection
     *   in barycentric coordinates
     * \return
     *   \c true if an intersection has been detected
     */
    FINLINE bool rayIntersect(const Point *positions, const Ray &ray, Float &u,
        Float &v, Float &t) const {
        return rayIntersect(
            positions[idx[0]], positions[idx[1]],
            positions[idx[2]], ray, u, v, t);
    }
#endif

};

#include <set>

class KDNode {
public:
    std::vector<KDTriangle*> triangles;
    KDNode* left;
    KDNode* right;
    AABB bbox;

    KDNode();
    KDNode(const KDNode& node);
    KDNode(std::vector<KDTriangle*>& tris);

    KDNode* build(std::vector<KDTriangle*>& tris, int depth) const;

    KDNode& operator= (const KDNode& rhs);

    int total();

    bool hit(const Ray& ray, float& t, float& tmin, hitInfo *info=NULL/*, ShadeRec& sr*/);
    int hitCount(const Ray& ray, float& t, float& tmin);
    bool isPointInside(Point p) ;
    //bool shadow_hit(KDNode* node, const Ray& ray, float& tmin, float Li) const;

private:
    inline bool aabbHit(const AABB& aabb, const Ray& ray) const {
        Float nearT, farT;
        return aabb.rayIntersect(ray, nearT, farT);
    }

};

} // namespace grain

#endif
