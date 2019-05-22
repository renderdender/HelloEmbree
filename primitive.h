#pragma once

#include <fstream>
#include <sstream>
#include <vector>

#include <embree3/rtcore.h>
#include <embree3/rtcore_ray.h>

#include "vector3.h"
#include "ray.h"
#include "bxdf.h"

/*************************************************    Declaration    **************************************************/

// TODO(?): Make a proper Transform class with matrices.
class Transform{

public:
    explicit Transform() : scale(), rotate(), translate() { }
    explicit Transform(Vec3 sc, Vec3 ro, Vec3 tr) : scale(sc), rotate(ro), translate(tr) { }

    Vec3 scale;
    Vec3 rotate;
    Vec3 translate;
};

class GeometricPrimitive{

public:
    GeometricPrimitive(Vec3 pos, Vec3 e, Vec3 c, BxDF_TYPE bxdf_type, Transform trans, bool if_inter = false, int prim_id = 0);

    virtual ~GeometricPrimitive() { };

    // Creating and committing the current object to Embree scene
    virtual int construct_embree_object(RTCDevice& rtc_device, RTCScene& rtc_scene) = 0;

    // Position before any transformation
    Vec3 position;
    Vec3 emission;
    Vec3 color;

    // TODO(?): Get rid of switch-case material system and replace it with proper OOP materials
    // BxDF type
    BxDF_TYPE bxdf;
    Transform transform;

    bool if_interpolated;
    int primitive_id;
};


class TriangleMesh final : public GeometricPrimitive{

public:
    explicit TriangleMesh(std::string file, Transform trans, Vec3 pos, Vec3 e, Vec3 c, BxDF_TYPE bxdf_type);

    ~TriangleMesh();

    // Creating and commiting the current object to Embree scene
    int construct_embree_object(RTCDevice& rtc_device, RTCScene& rtc_scene) override;

    std::string filename;

    std::vector<Vec3f> vertices;
    std::vector<Vec3f> normals;
    std::vector<int> indices_v;
    std::vector<int> indices_n;

    // Has to be freed in the destructor
    Vec3f* aligned_normals;
};


class Sphere final : public GeometricPrimitive{

public:
    explicit Sphere(Vec3::value_type rad, Transform trans, Vec3 pos, Vec3 e, Vec3 c, BxDF_TYPE bxdf_type);

    ~Sphere() = default;

    // User defined intersection functions for the Sphere primitive
    static void sphereBoundsFunc(const struct RTCBoundsFunctionArguments* args);
    static void sphereIntersectFunc(const RTCIntersectFunctionNArguments* args);
    static void sphereOccludedFunc(const RTCOccludedFunctionNArguments* args);

    // Creating and commiting the current object to Embree scene
    int construct_embree_object(RTCDevice& rtc_device, RTCScene& rtc_scene);

    Vec3::value_type radius;
};

/**************************************    Implementation: GeometricPrimitive    **************************************/

inline GeometricPrimitive::GeometricPrimitive(Vec3 pos, Vec3 e, Vec3 c, BxDF_TYPE bxdf_type, Transform trans, bool if_inter, int prim_id) :
    position(pos), emission(e), color(c), bxdf(bxdf_type), transform(trans), if_interpolated(if_inter), primitive_id(prim_id) { }

/****************************************    Implementation: TriangleMesh    ******************************************/

inline TriangleMesh::TriangleMesh(std::string file, Transform trans, Vec3 pos, Vec3 e, Vec3 c, BxDF_TYPE bxdf_type) :
    GeometricPrimitive(pos, e, c, bxdf_type, trans, true), filename(file) { }

inline TriangleMesh::~TriangleMesh(){
    delete [] aligned_normals;
}

// Construction of Embree object from the mesh
inline int TriangleMesh::construct_embree_object(RTCDevice& rtc_device, RTCScene& rtc_scene){
    // Reading .obj file
    std::ifstream in(filename, std::ios::in);
    if (!in) {
        std::cerr << "Cannot open " << filename << std::endl;
        return 1;

    }
    // Parsing .obj file.
    // It should contain:
    // - vertices (v)
    // - normals (vn)
    // - indices (f), formatted as idx//idx_normal idx//idx_normal idx//idx_normal
    // If your .obj files are formatted differently, rewrite the parser accordingly. Or just export your model with Blender.
    std::string line;
    while (std::getline(in, line)) {
        if (line.substr(0, 2) == "v ") {
            std::istringstream v(line.substr(2));
            float x, y, z;
            v >> x;
            v >> y;
            v >> z;
            Vec3f vertex = Vec3f(x, y, z) * Vec3f(this->transform.scale) + Vec3f(this->transform.translate);
            vertices.push_back(vertex);
        } else if (line.substr(0, 2) == "vn") {
            std::istringstream v(line.substr(2));
            float x, y, z;
            v >> x;
            v >> y;
            v >> z;
            Vec3f normal = Vec3f(x, y, z);
            normals.push_back(normal);
        } else if (line.substr(0, 2) == "f ") {
            int x, y, z;
            int x_n, y_n, z_n;
            const char *chh = line.c_str();
            sscanf(chh, "f %i//%i %i//%i %i//%i", &x, &x_n, &y, &y_n, &z, &z_n);

            indices_v.push_back(--x);
            indices_v.push_back(--y);
            indices_v.push_back(--z);
            indices_n.push_back(--x_n);
            indices_n.push_back(--y_n);
            indices_n.push_back(--z_n);
        }
    }

    /*
    // Manually calculated normals.
    indices_n = indices_v;
    normals = std::vector<Vec3f>(vertices.size());
    for (int i = 0; i < indices_v.size(); i += 3) {
        Vec3f edge_a = vertices[indices_v[i + 1]] - vertices[indices_v[i]];
        Vec3f edge_b = vertices[indices_v[i + 2]] - vertices[indices_v[i]];
        Vec3f current_normal = cross(edge_a, edge_b);

        normals[indices_v[i]] += current_normal;
        normals[indices_v[i + 1]] += current_normal;
        normals[indices_v[i + 2]] += current_normal;
    }

    for (auto &normal : normals) {
        normal = normalize(normal);
    }
    */

    unsigned int vertices_size = vertices.size();
    unsigned int indices_size = indices_v.size();

    // Initializing Embree geometry
    RTCGeometry mesh = rtcNewGeometry(rtc_device, RTC_GEOMETRY_TYPE_TRIANGLE);

    // Setting and filling the vertex buffer
    Vec3f *embree_vertices = (Vec3f *) rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
                                                               sizeof(Vec3f), vertices_size);
    for (int i = 0; i < vertices_size; ++i) {
        embree_vertices[i].x = vertices[i].x;
        embree_vertices[i].y = vertices[i].y;
        embree_vertices[i].z = vertices[i].z;
    }

    vertices.clear();
    aligned_normals = new Vec3f[vertices_size];

    // Setting and filling the index buffer
    Vec3i *triangles = (Vec3i *) rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
                                                         sizeof(Vec3i), indices_size / 3);
    for (int i = 0; i < indices_size; i += 3) {
        triangles[i / 3].x = indices_v[i];
        triangles[i / 3].y = indices_v[i + 1];
        triangles[i / 3].z = indices_v[i + 2];

        aligned_normals[indices_v[i]] = normals[indices_n[i]];
        aligned_normals[indices_v[i + 1]] = normals[indices_n[i + 1]];
        aligned_normals[indices_v[i + 2]] = normals[indices_n[i + 2]];
    }

    normals.clear();
    indices_n.clear();
    indices_v.clear();

    // Setting the buffer with normals (in order to get interpolated normals almost for free (since the uv-values are already calculated for intersections))
    rtcSetGeometryVertexAttributeCount(mesh, 1);

    rtcSetSharedGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3, &aligned_normals[0], 0,
                               sizeof(Vec3f), vertices_size);

    // Commiting the geometry
    rtcCommitGeometry(mesh);

    // Attaching it to the scene and getting the primitive's Id
    this->primitive_id = rtcAttachGeometry(rtc_scene, mesh);

    rtcReleaseGeometry(mesh);

    return 0;
}

/*******************************************    Implementation: Sphere    *********************************************/

inline Sphere::Sphere(Vec3::value_type rad, Transform trans, Vec3 pos, Vec3 e, Vec3 c, BxDF_TYPE bxdf_type) :
    GeometricPrimitive(pos, e, c, bxdf_type, trans, false), radius(rad) { }

// Bounding box construction routine
inline void Sphere::sphereBoundsFunc(const struct RTCBoundsFunctionArguments* args) {
    const Sphere* spheres = (const Sphere*) args->geometryUserPtr;
    RTCBounds* bounds_o = args->bounds_o;
    const Sphere& sphere = spheres[args->primID];

    bounds_o->lower_x = sphere.position.x-sphere.radius;
    bounds_o->lower_y = sphere.position.y-sphere.radius;
    bounds_o->lower_z = sphere.position.z-sphere.radius;
    bounds_o->upper_x = sphere.position.x+sphere.radius;
    bounds_o->upper_y = sphere.position.y+sphere.radius;
    bounds_o->upper_z = sphere.position.z+sphere.radius;

}

// Intersection routine
inline void Sphere::sphereIntersectFunc(const RTCIntersectFunctionNArguments* args) {
    int* valid = args->valid;
    void* ptr  = args->geometryUserPtr;
    RTCRayHit* rayhit = (RTCRayHit*)args->rayhit;
    unsigned int primID = args->primID;

    Ray t_ray;
    RTCRayHit_to_Ray(*rayhit, t_ray);

    assert(args->N == 1);

    const Sphere* spheres = (const Sphere*)ptr;
    const Sphere& sphere = spheres[primID];

    if (!valid[0]) return;

    const Vec3 op = sphere.position - t_ray.org;
    const double dop = dot(t_ray.dir, op);
    const double D = dop * dop - dot(op, op) + sphere.radius * sphere.radius;

    if (D < 0)
        return;

    const double sqrtD = std::sqrt(D);

    const double tmin = dop - sqrtD;
    if (t_ray.tnear < tmin && tmin < t_ray.tfar) {
        t_ray.tfar = tmin;

        rayhit->hit.u = 0.0f;
        rayhit->hit.v = 0.0f;
        rayhit->ray.tfar = tmin;
        rayhit->hit.geomID = sphere.primitive_id;
        rayhit->hit.primID = primID;
        Vec3 ng = t_ray.org + t_ray.dir * tmin - sphere.position;
        rayhit->hit.Ng_x = ng.x;
        rayhit->hit.Ng_y = ng.y;
        rayhit->hit.Ng_z = ng.z;
    }

    const double tmax = dop + sqrtD;
    if (t_ray.tnear < tmax && tmax < t_ray.tfar) {
        t_ray.tfar = tmax;

        rayhit->hit.u = 0.0f;
        rayhit->hit.v = 0.0f;
        rayhit->ray.tfar = tmax;
        rayhit->hit.geomID = sphere.primitive_id;
        rayhit->hit.primID = primID;
        Vec3 ng = t_ray.org + t_ray.dir * tmax - sphere.position;
        rayhit->hit.Ng_x = ng.x;
        rayhit->hit.Ng_y = ng.y;
        rayhit->hit.Ng_z = ng.z;
    }

    return;
}

// Occlusion routine
inline void Sphere::sphereOccludedFunc(const RTCOccludedFunctionNArguments* args) {
    // TODO: Implement the occlusion routine and make sure to use it for shadow rays or something similar
    // Hint: It's almost the same as sphereIntersectFunc
}

// Construction of Embree object from the analytically given sphere
inline int Sphere::construct_embree_object(RTCDevice& rtc_device, RTCScene& rtc_scene) {

    RTCGeometry geom = rtcNewGeometry(rtc_device, RTC_GEOMETRY_TYPE_USER);
    this->primitive_id = rtcAttachGeometry(rtc_scene, geom);

    rtcSetGeometryUserPrimitiveCount(geom, 1);
    rtcSetGeometryUserData(geom, this);

    rtcSetGeometryBoundsFunction(geom, sphereBoundsFunc, nullptr);
    rtcSetGeometryIntersectFunction(geom, sphereIntersectFunc);
    rtcSetGeometryOccludedFunction(geom, sphereOccludedFunc);
    rtcCommitGeometry(geom);
    rtcReleaseGeometry(geom);

    return 1;
}
