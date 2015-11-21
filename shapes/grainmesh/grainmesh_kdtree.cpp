#include "grainmesh_kdtree.h"

using namespace grain;

KDNode::KDNode() : triangles(std::vector<KDTriangle*>()), left(NULL), right(NULL), bbox(AABB())
{}

KDNode::KDNode(const KDNode& node) : triangles(node.triangles), left(node.left), right(node.right), bbox(node.bbox)
{}

KDNode::KDNode(std::vector<KDTriangle *> &tris) : triangles(std::vector<KDTriangle*>()), left(NULL), right(NULL)
{
    KDNode* node = build(tris, 0);
    triangles = node->triangles;
    left = node->left;
    right = node->right;
    bbox = node->bbox;
}

KDNode* KDNode::build(std::vector<KDTriangle *> &tris, int depth) const{
    KDNode* node = new KDNode();
    node->triangles = tris;
    node->left = NULL;
    node->right = NULL;
    node->bbox = AABB();

    if (tris.size() == 0)
        return node;


    if (tris.size() == 1) {
        //node->bbox = tris[0]->get_bounding_box();
        node->bbox = (tris.at(0))->getAABB();
        node->left = new KDNode();
        node->right = new KDNode();
        node->left->triangles = std::vector<KDTriangle*>();
        node->right->triangles = std::vector<KDTriangle*>();
        return node;
    }

    // get a bounding box surrounding all the triangles
    //node->bbox = tris[0]->get_bounding_box();
    node->bbox = (tris.at(0))->getAABB();

    for (unsigned int i=1; i < tris.size(); i++) {
        //node->bbox.expand(tris[i]->get_bounding_box());
        node->bbox.expandBy((tris.at(i))->getAABB());
    }

    Point midpt = Point(0,0,0);
    for (unsigned int i=0; i < tris.size(); i++) {
        // find midpoint of all triangles
        //midpt = midpt + (tris[i]->get_midpoint() * (1.0f / (float)tris.size()));
        midpt = midpt + tris[i]->get_midpoint();
    }
    midpt = midpt * (1.0f / (float)tris.size());

    std::vector<KDTriangle*> left_tris;
    std::vector<KDTriangle*> right_tris;
    //int axis = node->bbox.longest_axis();
    int axis = 0;
    if((node->bbox.max.y - node->bbox.min.y) > (node->bbox.max.x - node->bbox.min.x)
    && (node->bbox.max.y - node->bbox.min.y) > (node->bbox.max.z - node->bbox.min.z))
        axis = 1;
    if((node->bbox.max.z - node->bbox.min.z) > (node->bbox.max.x - node->bbox.min.x)
    && (node->bbox.max.z - node->bbox.min.z) > (node->bbox.max.y - node->bbox.min.y))
        axis = 2;

    for (unsigned int i=0; i < tris.size(); i++) {
        // split triangles based on their midpoints side of avg in longest axis
        switch (axis) {
            case 0:
                midpt.x >= tris[i]->get_midpoint().x ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
                break;
            case 1:
                midpt.y >= tris[i]->get_midpoint().y ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
                break;
            case 2:
                midpt.z >= tris[i]->get_midpoint().z ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
                break;
        }
    }

    if (left_tris.size() == 0 && right_tris.size() > 0) left_tris = right_tris;
    if (right_tris.size() == 0 && left_tris.size() > 0) right_tris = left_tris;

    // if 50% of triangles match, don't subdivide any more
    // actually when only 2 left in this node and on the midpt line
    int matches = 0;
    for (unsigned int i=0; i < left_tris.size(); i++) {
        for (unsigned int j=0; j < right_tris.size(); j++) {
            if (left_tris[i] == right_tris[j])
                matches++;
        }
    }

    if ((float)matches / left_tris.size() < 0.5 && (float)matches / right_tris.size() < 0.5) {
        // recurse down left and right sides
        node->left = build(left_tris, depth + 1);
        node->right = build(right_tris, depth + 1);
    }
    else {
        node->left = new KDNode();
        node->right = new KDNode();
        node->left->triangles = std::vector<KDTriangle*>();
        node->right->triangles = std::vector<KDTriangle*>();
    }

    return node;
}

KDNode& KDNode::operator=(const KDNode& rhs) {
    if (this == &rhs)
        return *this;

    bbox = rhs.bbox;
    left = rhs.left;
    right = rhs.right;
    for (unsigned int i=0; i < rhs.triangles.size(); i++) {
        triangles.push_back(rhs.triangles[i]);
    }

    return *this;
}


bool KDNode::hit(const Ray& ray, float& t, float& tmin, hitInfo* info/*, ShadeRec& sr*/) {
    //if (node->bbox.hit(ray)) {
    if(aabbHit(this->bbox, ray)) {
        //Normal normal;
        bool hit_tri = false;
        //Point hit_pt, local_hit_pt;
        bool hitl = false, hitr = false;

        if ((this->left != NULL && this->left->triangles.size() > 0) || (this->right != NULL && this->right->triangles.size() > 0)) {
            if(this->left != NULL && this->left->triangles.size() > 0)
                hitl = this->left->hit(ray, t, tmin, info/*, sr*/);
            if(this->right != NULL && this->right->triangles.size() > 0)
                hitr = this->right->hit(ray, t, tmin, info/*, sr*/);
            return hitl || hitr;
        }
        // leaf node
        else {
            float nearestT = 1.0/0.0;
            float to;
            KDTriangle *hitTriangle;
            for (size_t i=0; i < this->triangles.size(); i++) {
                if (this->triangles[i]->hit(ray, to, tmin/*, sr*/)) {
                    if (to<nearestT) {
                        hit_tri = true;
                        nearestT = to;
                        hitTriangle = this->triangles[i];
                    }
                    //ying sr.hit_obj = true;
                    //ying tmin = t;
                    //sr.ph = ray.o + t * ray.d;
                    //ying hit_pt = sr.ph;
                    //sr.nh = node->triangles[i]->get_normal();
                    //ying normal = sr.nh;
                    //ying sr.local_ph = node->triangles[i]->translate_to_local_coords(sr.ph);
                    //ying local_hit_pt = sr.local_ph;
                }
            }
            if (hit_tri) {
                /*
                sr.hit_obj = true;
                sr.t = tmin;
                sr.nh = normal;
                sr.ph = hit_pt;
                sr.local_ph = local_hit_pt;
                */

                if(nearestT < t) {
                    t = nearestT;
                    if(NULL != info) {
                        info->hitObject = (void*)hitTriangle;
                    }
                }
                return true;
            }
            return false;
        }
    }
    return false;
}

int KDNode::hitCount(const Ray& ray, float& t, float& tmin) {
    int count = 0;

    if(aabbHit(this->bbox, ray)) {

        if ((this->left != NULL && this->left->triangles.size() > 0) || (this->right != NULL && this->right->triangles.size() > 0)) {
            if(this->left != NULL && this->left->triangles.size() > 0)
                count += this->left->hitCount(ray, t, tmin);
            if(this->right != NULL && this->right->triangles.size() > 0)
                count += this->right->hitCount(ray, t, tmin);
            return count;
        }
        // leaf node
        else {
            float to;
            for (size_t i=0; i < this->triangles.size(); i++) {
                if (this->triangles[i]->hit(ray, to, tmin) && to>tmin) {
                    count++;
                }
            }
            return count;
        }
    }
    return count;
}

bool KDNode::isPointInside(Point p) {

    float tout = 1.0/0.0, tmin = 0.;
    return ((tmin = 0., tout = 1.0/0.0, hit(Ray(p, Vector(1,0,0), (Float)0), tout, tmin))
     && (tmin = 0., tout = 1.0/0.0, hit(Ray(p, Vector(-1,0,0), (Float)0), tout, tmin))
     && (tmin = 0., tout = 1.0/0.0, hit(Ray(p, Vector(0,1,0), (Float)0), tout, tmin))
     && (tmin = 0., tout = 1.0/0.0, hit(Ray(p, Vector(0,-1,0), (Float)0), tout, tmin))
     && (tmin = 0., tout = 1.0/0.0, hit(Ray(p, Vector(0,0,1), (Float)0), tout, tmin))
     && (tmin = 0., tout = 1.0/0.0, hit(Ray(p, Vector(0,0,-1), (Float)0), tout, tmin)));


    /*
    float t=1.0/0.0, tmin = 0.;
    int a = (hitCount(Ray(p, Vector(1,0,0), 0), t, tmin));
    int b = (hitCount(Ray(p, Vector(-1,0,0), 0), t, tmin));
    int c = (hitCount(Ray(p, Vector(0,1,0), 0), t, tmin));
    int d = (hitCount(Ray(p, Vector(0,-1,0), 0), t, tmin));
    int e = (hitCount(Ray(p, Vector(0,0,1), 0), t, tmin));
    int f = (hitCount(Ray(p, Vector(0,0,-1), 0), t, tmin));

    return (((hitCount(Ray(p, Vector(1,0,0), 0), t, tmin)) % 2 == 1) &&
           ((hitCount(Ray(p, Vector(-1,0,0), 0), t, tmin)) % 2 == 1)) &&
           (((hitCount(Ray(p, Vector(0,1,0), 0), t, tmin)) % 2 == 1) &&
           ((hitCount(Ray(p, Vector(0,-1,0), 0), t, tmin)) % 2 == 1)) &&
           (((hitCount(Ray(p, Vector(0,0,1), 0), t, tmin)) % 2 == 1) &&
           ((hitCount(Ray(p, Vector(0,0,-1), 0), t, tmin)) % 2 == 1));
           */
}

#if 0

bool KDNode::hit(KDNode* node, const Ray& ray, float& t, float& tmin, ShadeRec& sr) const {
    if (node->bbox.hit(ray)) {
        Normal normal;
        bool hit_tri = false;
        Point3D hit_pt, local_hit_pt;

        if (node->left != NULL && node->left->triangles.size() > 0 || node->right != NULL && node->right->triangles.size() > 0) {
            bool hitleft = hit(node->left, ray, t, tmin, sr);
            bool hitright = hit(node->right, ray, t, tmin, sr);
            return hitleft || hitright;
        }
        else {
            for (unsigned int i=0; i < node->triangles.size(); i++) {
                if (node->triangles[i]->hit(ray, t, tmin, sr)) {
                    hit_tri = true;
                    sr.hit_obj = true;
                    tmin = t;
                    //sr.ph = ray.o + t * ray.d;
                    hit_pt = sr.ph;
                    //sr.nh = node->triangles[i]->get_normal();
                    normal = sr.nh;
                    sr.local_ph = node->triangles[i]->translate_to_local_coords(sr.ph);
                    local_hit_pt = sr.local_ph;
                }
            }
            if (hit_tri) {
                sr.hit_obj = true;
                sr.t = tmin;
                sr.nh = normal;
                sr.ph = hit_pt;
                sr.local_ph = local_hit_pt;
                return true;
            }
            return false;
        }
    }
    return false;
}



bool KDNode::shadow_hit(KDNode* node, const Ray& ray, float& tmin, float Li) const {
    if (node->bbox.hit(ray)/* || bbox.inside(ray.o)*/) {
        float t;
        bool shadow_hit = false;

        if (node->left != NULL && node->left->triangles.size() > 0 || node->right != NULL && node->right->triangles.size() > 0) {
            bool hitleft = node->shadow_hit(node->left, ray, tmin, Li);
            bool hitright = node->shadow_hit(node->right, ray, tmin, Li);
            return hitleft || hitright;
        }
        else {
            for (unsigned int i=0; i < node->triangles.size(); i++) {
                if (node->triangles[i]->shadow_hit(ray, t, Li)) {
                    shadow_hit = true;
                    tmin = t;
                }
            }
            if (shadow_hit) {
                return true;
            }
            return false;
        }
    }
    return false;
}
#endif
