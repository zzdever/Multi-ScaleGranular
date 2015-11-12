
#if !defined(__GRAINMESH_H)
#define __GRAINMESH_H

#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

//class HairKDTree;

class GrainMesh : public Shape {
public:
    /// Construct a new HairShape instance given a properties object
    GrainMesh(const Properties &props);

    /// Unserialize from a binary data stream
    GrainMesh(Stream *stream, InstanceManager *manager);

    /// Serialize to a binary data stream
    void serialize(Stream *stream, InstanceManager *manager) const;

    // =============================================================
    //! @{ \name Access the internal vertex data
    // =============================================================

    /// Return the list of vertices underlying the hair shape
    const std::vector<Point> &getVertices() const;

    /**
     * Return a boolean list specifying whether a vertex
     * marks the beginning of a new fiber
     */
    const std::vector<bool> &getStartFiber() const;

    //! @}
    // =============================================================

    // =============================================================
    //! @{ \name Implementation of the \ref Shape interface
    // =============================================================

    bool rayIntersect(const Ray &ray, Float mint,
            Float maxt, Float &t, void *temp) const;

    bool rayIntersect(const Ray &ray, Float mint, Float maxt) const;

    void fillIntersectionRecord(const Ray &ray,
        const void *temp, Intersection &its) const;

    ref<TriMesh> createTriMesh();

    const KDTreeBase<AABB> *getKDTree() const;

    AABB getAABB() const;

    Float getSurfaceArea() const;

    size_t getPrimitiveCount() const;

    size_t getEffectivePrimitiveCount() const;

    //! @}
    // =============================================================

    /// Return a human-readable representation
    std::string toString() const;

    MTS_DECLARE_CLASS()
private:
    //ref<HairKDTree> m_kdtree;

protected:
    std::string m_name;
};

MTS_NAMESPACE_END

#endif /* __GRAINMESH_H */
