
#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/fresolver.h>

#include "../../shapes/grainmesh/tempinfo.h"

#include <fstream>


MTS_NAMESPACE_BEGIN

#define OUTPUT_CONDITION (control_timer->getMicroseconds()%1000 == 0)


static StatsCounter avgPathLength("Path tracer", "Average path length", EAverage);

static StatsCounter count("Path tracer", "marker", ENumberValue);


class MultiScale : public MonteCarloIntegrator {
public:
    MultiScale(const Properties &props)
        : MonteCarloIntegrator(props) { }

    /// Unserialize from a binary data stream
    MultiScale(Stream *stream, InstanceManager *manager)
        : MonteCarloIntegrator(stream, manager) { }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        bool scattered = false;

        MediumSamplingRecord mRec;
        bool nullChain = false;


        MODES_T rendermode;
        // backup the address, its.temp will be modified in rayIntersect()
        ScatterometerReturnInfo* scatterometerReturnInfo = static_cast<ScatterometerReturnInfo*>(its.temp);
        bool isScatterometer = false;
        bool isIntersected = false;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);
        ray.mint = Epsilon;

        Spectrum throughput(1.0f);
        Float eta = 1.0f;

        if (m_maxDepth == 1)
            rRec.type &= RadianceQueryRecord::EEmittedRadiance;


        GrainIntersectionInfo* info = static_cast<GrainIntersectionInfo*>(its.temp);

        isIntersected = its.isValid();
        if(its.isValid()) {
            // rendermode can be abtained after at least one valid intersection
            rendermode = info->rendermode;

            if(SCATTEROMETER==rendermode) {
                new(scatterometerReturnInfo) ScatterometerReturnInfo;

                scatterometerReturnInfo->intersection = its.p;
                scatterometerReturnInfo->outDirection = ray.d; // record the out direction for the last valid intersection

                /* mark as scatterometer */
                isScatterometer = true;
            }
        }


        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {

#define OUT 0
#if OUT
            std::cout<<"[multiscale1]: "<<rRec.toString()<<std::endl;
#endif

            if (!its.isValid()) {
                /* If no intersection could be found, potentially return
                   radiance from a environment luminaire if it exists */
                if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered))
                    Li += throughput * scene->evalEnvironment(ray);

#if OUT
                    std::cout<<"break out 0"<<std::endl;
#endif

                break;
            } else {
                isIntersected = true;
            }


            if(VPT==rendermode)
            {

#if OUT
                std::cout<<"!!!!!!!!!! in VPT!!!!!!!"<<std::endl;
#endif


                /* ==================================================================== */
                /*                 Radiative Transfer Equation sampling                 */
                /* ==================================================================== */
                if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
                    /* Sample the integral
                       \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
                    */


#if OUT
                    std::cout<<"!!!!!!!!!! in medium!!!!!!!"<<std::endl;
#endif


                    const PhaseFunction *phase = rRec.medium->getPhaseFunction();

                    throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

                    /* ==================================================================== */
                    /*                     Direct illumination sampling                     */
                    /* ==================================================================== */

                    /* Estimate the single scattering component if this is requested */
                    if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
                        DirectSamplingRecord dRec(mRec.p, mRec.time);
                        int maxInteractions = m_maxDepth - rRec.depth - 1;

                        Spectrum value = scene->sampleAttenuatedEmitterDirect(
                                dRec, rRec.medium, maxInteractions,
                                rRec.nextSample2D(), rRec.sampler);

                        if (!value.isZero())
                            Li += throughput * value * phase->eval(
                                    PhaseFunctionSamplingRecord(mRec, -ray.d, dRec.d));
                    }


                    /* Stop if multiple scattering was not requested, or if the path gets too long */
                    if ((rRec.depth + 1 >= m_maxDepth && m_maxDepth > 0) ||
                        !(rRec.type & RadianceQueryRecord::EIndirectMediumRadiance))
                        break;

                    /* ==================================================================== */
                    /*             Phase function sampling / Multiple scattering            */
                    /* ==================================================================== */

                    PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                    Float phaseVal = phase->sample(pRec, rRec.sampler);
                    if (phaseVal == 0)
                        break;
                    throughput *= phaseVal;

                    /* Trace a ray in this direction */
                    ray = Ray(mRec.p, pRec.wo, ray.time);
                    ray.mint = 0;
                    scene->rayIntersect(ray, its);
                    nullChain = false;
                    scattered = true;
                } else {
                    /* Sample
                        tau(x, y) * (Surface integral). This happens with probability mRec.pdfFailure
                        Account for this and multiply by the proper per-color-channel transmittance.
                    */

                    if (rRec.medium)
                        throughput *= mRec.transmittance / mRec.pdfFailure;

                    if (!its.isValid()) {
                        /* If no intersection could be found, possibly return
                           attenuated radiance from a background luminaire */
                        if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                            && (!m_hideEmitters || scattered)) {
                            Spectrum value = throughput * scene->evalEnvironment(ray);
                            if (rRec.medium)
                                value *= rRec.medium->evalTransmittance(ray);
                            Li += value;
                        }
                        break;
                    }

                    /* Possibly include emitted radiance if requested */
                    if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                        && (!m_hideEmitters || scattered))
                        Li += throughput * its.Le(-ray.d);

                    /* Include radiance from a subsurface integrator if requested */
                    if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                        Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

                    /* Prevent light leaks due to the use of shading normals */
                    Float wiDotGeoN = -dot(its.geoFrame.n, ray.d),
                          wiDotShN  = Frame::cosTheta(its.wi);
                    if (m_strictNormals && wiDotGeoN * wiDotShN < 0)
                        break;

                    /* ==================================================================== */
                    /*                     Direct illumination sampling                     */
                    /* ==================================================================== */

                    const BSDF *bsdf = its.getBSDF(ray);

                    /* Estimate the direct illumination if this is requested */
                    if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                            (bsdf->getType() & BSDF::ESmooth)) {
                        DirectSamplingRecord dRec(its);
                        int maxInteractions = m_maxDepth - rRec.depth - 1;

                        Spectrum value = scene->sampleAttenuatedEmitterDirect(
                                dRec, its, rRec.medium, maxInteractions,
                                rRec.nextSample2D(), rRec.sampler);

                        if (!value.isZero()) {
                            /* Allocate a record for querying the BSDF */
                            BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
                            bRec.sampler = rRec.sampler;

                            Float woDotGeoN = dot(its.geoFrame.n, dRec.d);
                            /* Prevent light leaks due to the use of shading normals */
                            if (!m_strictNormals ||
                                woDotGeoN * Frame::cosTheta(bRec.wo) > 0)
                                Li += throughput * value * bsdf->eval(bRec);
                        }
                    }

                    /* ==================================================================== */
                    /*                   BSDF sampling / Multiple scattering                */
                    /* ==================================================================== */

                    /* Sample BSDF * cos(theta) */
                    BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
                    Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
                    if (bsdfVal.isZero())
                        break;

                    /* Recursively gather indirect illumination? */
                    int recursiveType = 0;
                    if ((rRec.depth + 1 < m_maxDepth || m_maxDepth < 0) &&
                        (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                        recursiveType |= RadianceQueryRecord::ERadianceNoEmission;

                    /* Recursively gather direct illumination? This is a bit more
                       complicated by the fact that this integrator can create connection
                       through index-matched medium transitions (ENull scattering events) */
                    if ((rRec.depth < m_maxDepth || m_maxDepth < 0) &&
                        (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) &&
                        (bRec.sampledType & BSDF::EDelta) &&
                        (!(bRec.sampledType & BSDF::ENull) || nullChain)) {
                        recursiveType |= RadianceQueryRecord::EEmittedRadiance;
                        nullChain = true;
                    } else {
                        nullChain &= bRec.sampledType == BSDF::ENull;
                    }

                    /* Potentially stop the recursion if there is nothing more to do */
                    if (recursiveType == 0)
                        break;
                    rRec.type = recursiveType;

                    /* Prevent light leaks due to the use of shading normals */
                    const Vector wo = its.toWorld(bRec.wo);
                    Float woDotGeoN = dot(its.geoFrame.n, wo);
                    if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
                        break;

                    /* Keep track of the throughput, medium, and relative
                       refractive index along the path */
                    throughput *= bsdfVal;
                    eta *= bRec.eta;
                    if (its.isMediumTransition())
                        rRec.medium = its.getTargetMedium(wo);

                    /* In the next iteration, trace a ray in this direction */
                    ray = Ray(its.p, wo, ray.time);
                    scene->rayIntersect(ray, its);
                    scattered |= bRec.sampledType != BSDF::ENull;
                }

#if 0
                // VPT
                ray.o = ray.o + its.t*ray.d;

                double theta = (1.0*random()/RAND_MAX - 0.5) * M_PI; // -pi/2 ~ pi/2
                double phi = (1.0*random()/RAND_MAX - 0.5) * M_PI; // -pi/2 ~ pi/2
                ray.d = normalize(Vector(cos(phi)*cos(theta), sin(phi), -cos(phi)*sin(theta)));
                Point out = ray.o + ray.d;

                /* Trace a ray in this direction */
                ray = Ray(out, ray.d, ray.time);

                ray.mint = 0;
                scene->rayIntersect(ray, its);
                //nullChain = false;
                scattered = true;
#endif
            } // ~VPT
            else // EPT, SCATTEROMETER
            {
                const BSDF *bsdf = its.getBSDF(ray);

                /* Possibly include emitted radiance if requested */
                if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered))
                    Li += throughput * its.Le(-ray.d);

                /* Include radiance from a subsurface scattering model if requested */
                if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                    Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

                if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
                    || (m_strictNormals && dot(ray.d, its.geoFrame.n)
                        * Frame::cosTheta(its.wi) >= 0)) {

                    /* Only continue if:
                       1. The current path length is below the specifed maximum
                       2. If 'strictNormals'=true, when the geometric and shading
                          normals classify the incident direction to the same side */
#if OUT
                    std::cout<<"break out 1"<<std::endl;
#endif
                    break;
                }

                /* ==================================================================== */
                /*                     Direct illumination sampling                     */
                /* ==================================================================== */

                /* Estimate the direct illumination if this is requested */
                DirectSamplingRecord dRec(its);

                if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                    (bsdf->getType() & BSDF::ESmooth)) {
                    Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
                    if (!value.isZero()) {
                        const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                        /* Allocate a record for querying the BSDF */
                        BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

                        /* Evaluate BSDF * cos(theta) */
                        const Spectrum bsdfVal = bsdf->eval(bRec);

                        /* Prevent light leaks due to the use of shading normals */
                        if (!bsdfVal.isZero() && (!m_strictNormals
                                || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

                            /* Calculate prob. of having generated that direction
                               using BSDF sampling */
                            Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                                ? bsdf->pdf(bRec) : 0;

                            /* Weight using the power heuristic */
                            Float weight = miWeight(dRec.pdf, bsdfPdf);
                            Li += throughput * value * bsdfVal * weight;
                        }
                    }
                }


#if OUT
                    std::cout<<"before bsdf"<<std::endl;
#endif


                /* ==================================================================== */
                /*                            BSDF sampling                             */
                /* ==================================================================== */

                /* Sample BSDF * cos(theta) */
                Float bsdfPdf;
                BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
                Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
                if (bsdfWeight.isZero())
                    break;

#if OUT
                    std::cout<<"break out 2"<<std::endl;
#endif

                scattered |= bRec.sampledType != BSDF::ENull;

                /* Prevent light leaks due to the use of shading normals */
                const Vector wo = its.toWorld(bRec.wo);
                Float woDotGeoN = dot(its.geoFrame.n, wo);
                if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                    break;

#if OUT
                    std::cout<<"break out 3"<<std::endl;
                    std::cout<<"next ray direction: "<<wo.toString()<<std::endl;
#endif

                bool hitEmitter = false;
                Spectrum value;

                /* Record this out direction for a valid intersection */
                if(isScatterometer) {
                    scatterometerReturnInfo->intersection = its.p;
                    scatterometerReturnInfo->outDirection = wo;
                }

                /* Trace a ray in this direction */
                ray = Ray(its.p, wo, ray.time);
                if (scene->rayIntersect(ray, its)) {
                    /* Intersected something - check if it was a luminaire */
                    if (its.isEmitter()) {
                        value = its.Le(-ray.d);
                        dRec.setQuery(ray, its);
                        hitEmitter = true;
                    }
                } else {
                    /* Intersected nothing -- perhaps there is an environment map? */
                    const Emitter *env = scene->getEnvironmentEmitter();

                    if (env) {
                        if (m_hideEmitters && !scattered)
                            break;

                        value = env->evalEnvironment(ray);
                        if (!env->fillDirectSamplingRecord(dRec, ray))
                            break;
                        hitEmitter = true;
                    } else {
                        break;
                    }
                }

                /* Keep track of the throughput and relative
                   refractive index along the path */
                throughput *= bsdfWeight;
                eta *= bRec.eta;

                /* If a luminaire was hit, estimate the local illumination and
                   weight using the power heuristic */
                if (hitEmitter &&
                    (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                    /* Compute the prob. of generating that direction using the
                       implemented direct illumination sampling technique */
                    const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                        scene->pdfEmitterDirect(dRec) : 0;
                    Li += throughput * value * miWeight(bsdfPdf, lumPdf);
                }

                /* ==================================================================== */
                /*                         Indirect illumination                        */
                /* ==================================================================== */

                /* Set the recursive query type. Stop if no surface was hit by the
                   BSDF sample or if indirect illumination was not requested */
                if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                    break;
                rRec.type = RadianceQueryRecord::ERadianceNoEmission;

            } // ~EPT


            if (rRec.depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (rRec.nextSample1D() >= q)
                    break;
                throughput /= q;
            }
        } // while




        if(isScatterometer) {
#if OUT
            std::cout<<"[multiscale3]: "<<ray.toString()<<std::endl;

            std::cout<<"[multiscale-returninfo]: size: "<<sizeof(ScatterometerReturnInfo)
                    <<" "<<scatterometerReturnInfo->intersection.toString()<<std::endl<<scatterometerReturnInfo->outDirection.toString()<<std::endl;

#endif
            //rRec.its.t = 1.0f; // set it to valid
            //rRec.its.wi = outDir; // the out direction for last valid interseciton

            // do not store statistics
            return Li;
        }

        /*
        if(isScatterometer && isIntersected) {
            rRec.its.t = 1.0f; // set it to valid
            rRec.its.wi = outDir; // the out direction for last valid interseciton
        }
        */




        /* Store statistics */
        avgPathLength.incrementBase();
        avgPathLength += rRec.depth;

        return Li;
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MonteCarloIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MultiScale[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};






#if 0
class MultiScale : public MonteCarloIntegrator {
public:
    MultiScale(const Properties &props) : MonteCarloIntegrator(props) { }

    /// Unserialize from a binary data stream
    MultiScale(Stream *stream, InstanceManager *manager)
     : MonteCarloIntegrator(stream, manager) { }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        MediumSamplingRecord mRec;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        bool nullChain = true, scattered = false;
        Float eta = 1.0f;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);
        Spectrum throughput(1.0f);

        if (m_maxDepth == 1)
            rRec.type &= RadianceQueryRecord::EEmittedRadiance;

        /**
         * Note: the logic regarding maximum path depth may appear a bit
         * strange. This is necessary to get this integrator's output to
         * exactly match the output of other integrators under all settings
         * of this parameter.
         */

        ref<Timer> control_timer = new Timer();
        bool mm=false;
        std::vector<Point> trajectory;
        std::vector<Point> spheres;
        std::vector<int> spheres_indeces;
        std::vector<Point> origins;


        std::string logss;
        logss+="\n";
        count+=1;
        int iid = count.getValue();

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            if (!its.isValid()) {
                /* If no intersection could be found, possibly return
                   attenuated radiance from a background luminaire */
                if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered)) {
                    Spectrum value = throughput * scene->evalEnvironment(ray);
                    if (rRec.medium)
                        value *= rRec.medium->evalTransmittance(ray);
                    Li += value;
                }
                break;
            }

            /* Possibly include emitted radiance if requested */
            if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                && (!m_hideEmitters || scattered))
                Li += throughput * its.Le(-ray.d);


            /* ==================================================================== */
            /*                     Direct illumination sampling                     */
            /* ==================================================================== */

            /* Estimate the single scattering component if this is requested */
            if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
                DirectSamplingRecord dRec(mRec.p, mRec.time);
                int maxInteractions = m_maxDepth - rRec.depth - 1;

                Spectrum value = scene->sampleAttenuatedEmitterDirect(
                        dRec, rRec.medium, maxInteractions,
                        rRec.nextSample2D(), rRec.sampler);

                if (!value.isZero())
                    Li += throughput * value * 0;
                            //* phase->eval(
                            //PhaseFunctionSamplingRecord(mRec, -ray.d, dRec.d));
            }


            /* Stop if the path gets too long */
            if ((rRec.depth + 1 >= m_maxDepth && m_maxDepth > 0))
                break;

            /* ==================================================================== */
            /*             Phase function sampling / Multiple scattering            */
            /* ==================================================================== */

            /*
            PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
            Float phaseVal = phase->sample(pRec, rRec.sampler);
            if (phaseVal == 0)
                break;
            throughput *= phaseVal;
            */



            Li += throughput * its.color * 0.1;

            // EPT
            GrainIntersectionInfo* info = static_cast<GrainIntersectionInfo*>(its.temp);

#if 0
            {
                gotint.push_back(Point(info->ix, info->iy, info->iz));
                origins.push_back(Point(info->ox, info->oy, info->oz));
                tts.push_back(info->t);
                tts2.push_back(rRec.its.t);
                spheres.push_back(Point(info->x, info->y, info->z));
                spheres_indeces.push_back(info->name);
                trajectory.push_back(its.p);

                if(0 & OUTPUT_CONDITION)
                    Log(EWarn, "temppppp:%f, %f, %f, %f", info->center.x, info->center.y, info->center.z, info->radius);
            }
#endif


#if 0
            // VPT
            ray.o = ray.o + its.t*ray.d;

            double theta = (1.0*random()/RAND_MAX - 0.5) * M_PI; // -pi/2 ~ pi/2
            double phi = (1.0*random()/RAND_MAX - 0.5) * M_PI; // -pi/2 ~ pi/2
            ray.d = normalize(Vector(cos(phi)*cos(theta), sin(phi), -cos(phi)*sin(theta)));
            Point out = ray.o + ray.d;

            /* Trace a ray in this direction */
            ray = Ray(out, ray.d, ray.time);

            ray.mint = 0;
            scene->rayIntersect(ray, its);
            nullChain = false;
            scattered = true;

#endif


            Point intersectionPoint = ray.o + rRec.its.t * ray.d;
            Point out = ray.o +
                    (rRec.its.t + 2 * (dot(ray.d,
                                           (Vector(info->center.x, info->center.y, info->center.z)
                                            - Vector(intersectionPoint.x, intersectionPoint.y, intersectionPoint.z)))))
                    * ray.d;


            Normal n = its.geoFrame.n;
            ray.d = -2 * dot(ray.d, normalize(n)) * normalize(n) + ray.d;

            /* Trace a ray in this direction */
            ray = Ray(out, ray.d, ray.time);
            //ray = Ray(mRec.p, pRec.wo, ray.time);
            ray.mint = 0;
            scene->rayIntersect(ray, its);
            nullChain = false;
            scattered = true;






#if 0

            /* ==================================================================== */
            /*                 Radiative Transfer Equation sampling                 */
            /* ==================================================================== */

            if (rRec.medium && rRec.medium->sampleDistance(Ray(ray, 0, its.t), mRec, rRec.sampler)) {
                /* Sample the integral
                   \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
                */
                const PhaseFunction *phase = rRec.medium->getPhaseFunction();

                throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

                /* ==================================================================== */
                /*                     Direct illumination sampling                     */
                /* ==================================================================== */

                /* Estimate the single scattering component if this is requested */
                if (rRec.type & RadianceQueryRecord::EDirectMediumRadiance) {
                    DirectSamplingRecord dRec(mRec.p, mRec.time);
                    int maxInteractions = m_maxDepth - rRec.depth - 1;

                    Spectrum value = scene->sampleAttenuatedEmitterDirect(
                            dRec, rRec.medium, maxInteractions,
                            rRec.nextSample2D(), rRec.sampler);

                    if (!value.isZero())
                        Li += throughput * value * phase->eval(
                                PhaseFunctionSamplingRecord(mRec, -ray.d, dRec.d));
                }

                /* Stop if multiple scattering was not requested, or if the path gets too long */
                if ((rRec.depth + 1 >= m_maxDepth && m_maxDepth > 0) ||
                    !(rRec.type & RadianceQueryRecord::EIndirectMediumRadiance))
                    break;

                /* ==================================================================== */
                /*             Phase function sampling / Multiple scattering            */
                /* ==================================================================== */

                PhaseFunctionSamplingRecord pRec(mRec, -ray.d);
                Float phaseVal = phase->sample(pRec, rRec.sampler);
                if (phaseVal == 0)
                    break;
                throughput *= phaseVal;

                /* Trace a ray in this direction */
                ray = Ray(mRec.p, pRec.wo, ray.time);
                ray.mint = 0;
                scene->rayIntersect(ray, its);
                nullChain = false;
                scattered = true;
            } else {
                /* Sample
                    tau(x, y) * (Surface integral). This happens with probability mRec.pdfFailure
                    Account for this and multiply by the proper per-color-channel transmittance.
                */

                if (rRec.medium)
                    throughput *= mRec.transmittance / mRec.pdfFailure;

                if (!its.isValid()) {
                    /* If no intersection could be found, possibly return
                       attenuated radiance from a background luminaire */
                    if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                        && (!m_hideEmitters || scattered)) {
                        Spectrum value = throughput * scene->evalEnvironment(ray);
                        if (rRec.medium)
                            value *= rRec.medium->evalTransmittance(ray);
                        Li += value;
                    }
                    break;
                }


                /* Possibly include emitted radiance if requested */
                if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered))
                    Li += throughput * its.Le(-ray.d);

                /* Include radiance from a subsurface integrator if requested */
                if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                    Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

                /* Prevent light leaks due to the use of shading normals */
                Float wiDotGeoN = -dot(its.geoFrame.n, ray.d),
                      wiDotShN  = Frame::cosTheta(its.wi);
                if (m_strictNormals && wiDotGeoN * wiDotShN < 0)
                    break;


                /* ==================================================================== */
                /*                     Direct illumination sampling                     */
                /* ==================================================================== */

                const BSDF *bsdf = its.getBSDF(ray);

                /* Estimate the direct illumination if this is requested */
                if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                        (bsdf->getType() & BSDF::ESmooth)) {
                    DirectSamplingRecord dRec(its);
                    int maxInteractions = m_maxDepth - rRec.depth - 1;

                    Spectrum value = scene->sampleAttenuatedEmitterDirect(
                            dRec, its, rRec.medium, maxInteractions,
                            rRec.nextSample2D(), rRec.sampler);

                    if (!value.isZero()) {
                        /* Allocate a record for querying the BSDF */
                        BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));
                        bRec.sampler = rRec.sampler;

                        Float woDotGeoN = dot(its.geoFrame.n, dRec.d);
                        /* Prevent light leaks due to the use of shading normals */
                        if (!m_strictNormals ||
                            woDotGeoN * Frame::cosTheta(bRec.wo) > 0)
                            Li += throughput * value * bsdf->eval(bRec);
                    }
                }

                /* ==================================================================== */
                /*                   BSDF sampling / Multiple scattering                */
                /* ==================================================================== */

                /* Sample BSDF * cos(theta) */
                BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
                Spectrum bsdfVal = bsdf->sample(bRec, rRec.nextSample2D());
                if (bsdfVal.isZero())
                    break;

                /* Recursively gather indirect illumination? */
                int recursiveType = 0;
                if ((rRec.depth + 1 < m_maxDepth || m_maxDepth < 0) &&
                    (rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                    recursiveType |= RadianceQueryRecord::ERadianceNoEmission;

                /* Recursively gather direct illumination? This is a bit more
                   complicated by the fact that this integrator can create connection
                   through index-matched medium transitions (ENull scattering events) */
                if ((rRec.depth < m_maxDepth || m_maxDepth < 0) &&
                    (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance) &&
                    (bRec.sampledType & BSDF::EDelta) &&
                    (!(bRec.sampledType & BSDF::ENull) || nullChain)) {
                    recursiveType |= RadianceQueryRecord::EEmittedRadiance;
                    nullChain = true;
                } else {
                    nullChain &= bRec.sampledType == BSDF::ENull;
                }

                /* Potentially stop the recursion if there is nothing more to do */
                if (recursiveType == 0)
                    break;
                rRec.type = recursiveType;

                /* Prevent light leaks due to the use of shading normals */
                const Vector wo = its.toWorld(bRec.wo);
                Float woDotGeoN = dot(its.geoFrame.n, wo);
                if (woDotGeoN * Frame::cosTheta(bRec.wo) <= 0 && m_strictNormals)
                    break;

                /* Keep track of the throughput, medium, and relative
                   refractive index along the path */
                throughput *= bsdfVal;
                eta *= bRec.eta;
                if (its.isMediumTransition())
                    rRec.medium = its.getTargetMedium(wo);

                /* In the next iteration, trace a ray in this direction */
                ray = Ray(its.p, wo, ray.time);
                scene->rayIntersect(ray, its);
                scattered |= bRec.sampledType != BSDF::ENull;
            }

#endif

            if (rRec.depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (rRec.nextSample1D() >= q)
                    break;
                throughput /= q;
            }
        } // while

        avgPathLength.incrementBase();
        avgPathLength += rRec.depth;



#if 0
        if(0 && iid%1000==0 & trajectory.size()>1){
            std::stringstream ss;
            ss<<"\n";
            ss<<iid<<"\n";
            ss<<logss;
            ss<<"m_maxDepth: "<<m_maxDepth<<", m_rrDepth: "<<m_rrDepth<< ", rRec.depth: "<<rRec.depth<<"\n";
            ss<<"=============\n";
            Point p;
            for(int i=0; i<trajectory.size(); i++) {
                p = trajectory.at(i);
                ss<<p.x<<", "<< p.y<<", "<< p.z<<" || ";
            }
            ss<<"=============\n";
            Log(EWarn, ss.str().c_str());
        }
#endif

        return Li;
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MonteCarloIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MultiScale[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};

#endif




#if 0
static StatsCounter avgPathLength("Path tracer", "Average path length", EAverage);

class MultiScale : public MonteCarloIntegrator {
public:
    MultiScale(const Properties &props)
        : MonteCarloIntegrator(props) { }

    /// Unserialize from a binary data stream
    MultiScale(Stream *stream, InstanceManager *manager)
        : MonteCarloIntegrator(stream, manager) { }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        /* Some aliases and local variables */
        const Scene *scene = rRec.scene;
        Intersection &its = rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);

        bool scattered = false;

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        rRec.rayIntersect(ray);
        ray.mint = Epsilon;

        Spectrum throughput(1.0f);
        Float eta = 1.0f;

        while (rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            if (!its.isValid()) {
                /* If no intersection could be found, potentially return
                   radiance from a environment luminaire if it exists */
                if ((rRec.type & RadianceQueryRecord::EEmittedRadiance)
                    && (!m_hideEmitters || scattered))
                    Li += throughput * scene->evalEnvironment(ray);
                break;
            }

            const BSDF *bsdf = its.getBSDF(ray);

            /* Possibly include emitted radiance if requested */
            if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance)
                && (!m_hideEmitters || scattered))
                Li += throughput * its.Le(-ray.d);

            /* Include radiance from a subsurface scattering model if requested */
            if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                Li += throughput * its.LoSub(scene, rRec.sampler, -ray.d, rRec.depth);

            if ((rRec.depth >= m_maxDepth && m_maxDepth > 0)
                || (m_strictNormals && dot(ray.d, its.geoFrame.n)
                    * Frame::cosTheta(its.wi) >= 0)) {

                /* Only continue if:
                   1. The current path length is below the specifed maximum
                   2. If 'strictNormals'=true, when the geometric and shading
                      normals classify the incident direction to the same side */
                break;
            }



            Li = throughput*its.color;



#if 0
            /* ==================================================================== */
            /*                     Direct illumination sampling                     */
            /* ==================================================================== */

            /* Estimate the direct illumination if this is requested */
            DirectSamplingRecord dRec(its);

            if (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance &&
                (bsdf->getType() & BSDF::ESmooth)) {
                Spectrum value = scene->sampleEmitterDirect(dRec, rRec.nextSample2D());
                if (!value.isZero()) {
                    const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

                    /* Allocate a record for querying the BSDF */
                    BSDFSamplingRecord bRec(its, its.toLocal(dRec.d), ERadiance);

                    /* Evaluate BSDF * cos(theta) */
                    const Spectrum bsdfVal = bsdf->eval(bRec);

                    /* Prevent light leaks due to the use of shading normals */
                    if (!bsdfVal.isZero() && (!m_strictNormals
                            || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0)) {

                        /* Calculate prob. of having generated that direction
                           using BSDF sampling */
                        Float bsdfPdf = (emitter->isOnSurface() && dRec.measure == ESolidAngle)
                            ? bsdf->pdf(bRec) : 0;

                        /* Weight using the power heuristic */
                        Float weight = miWeight(dRec.pdf, bsdfPdf);
                        Li += throughput * value * bsdfVal * weight;
                    }
                }
            }
#endif
#if 0
            /* ==================================================================== */
            /*                            BSDF sampling                             */
            /* ==================================================================== */

            /* Sample BSDF * cos(theta) */
            Float bsdfPdf;
            BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
            Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
            if (bsdfWeight.isZero())
                break;

            scattered |= bRec.sampledType != BSDF::ENull;

            /* Prevent light leaks due to the use of shading normals */
            const Vector wo = its.toWorld(bRec.wo);
            Float woDotGeoN = dot(its.geoFrame.n, wo);
            if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                break;

            bool hitEmitter = false;
            Spectrum value;

            /* Trace a ray in this direction */
            ray = Ray(its.p, wo, ray.time);
            if (scene->rayIntersect(ray, its)) {
                /* Intersected something - check if it was a luminaire */
                if (its.isEmitter()) {
                    value = its.Le(-ray.d);
                    dRec.setQuery(ray, its);
                    hitEmitter = true;
                }
            } else {
                /* Intersected nothing -- perhaps there is an environment map? */
                const Emitter *env = scene->getEnvironmentEmitter();

                if (env) {
                    if (m_hideEmitters && !scattered)
                        break;

                    value = env->evalEnvironment(ray);
                    if (!env->fillDirectSamplingRecord(dRec, ray))
                        break;
                    hitEmitter = true;
                } else {
                    break;
                }
            }

            /* Keep track of the throughput and relative
               refractive index along the path */
            throughput *= bsdfWeight;
            eta *= bRec.eta;

            /* If a luminaire was hit, estimate the local illumination and
               weight using the power heuristic */
            if (hitEmitter &&
                (rRec.type & RadianceQueryRecord::EDirectSurfaceRadiance)) {
                /* Compute the prob. of generating that direction using the
                   implemented direct illumination sampling technique */
                const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ?
                    scene->pdfEmitterDirect(dRec) : 0;
                Li += throughput * value * miWeight(bsdfPdf, lumPdf);
            }
#endif
            /* ==================================================================== */
            /*                         Indirect illumination                        */
            /* ==================================================================== */

            /* Set the recursive query type. Stop if no surface was hit by the
               BSDF sample or if indirect illumination was not requested */
            if (!its.isValid() || !(rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                break;
            rRec.type = RadianceQueryRecord::ERadianceNoEmission;

            if (rRec.depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (rRec.nextSample1D() >= q)
                    break;
                throughput /= q;
            }

        }

        /* Store statistics */
        avgPathLength.incrementBase();
        avgPathLength += rRec.depth;

        return Li;
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MonteCarloIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MultiScale[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
};

#endif

MTS_IMPLEMENT_CLASS_S(MultiScale, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(MultiScale, "Multiscale rendering of granular materials");
MTS_NAMESPACE_END

