#pragma once

#include <embree3/rtcore_ray.h>

#include "vector3.h"

/*************************************************    Declaration    **************************************************/

class Ray {

public:
    // Initialize with default values
    explicit Ray();
    // Initialize with 2 vectors
    explicit Ray(Vec3 o, Vec3 d);

    ~Ray() = default;

    // Default copy constructor and assignment operator
    Ray(const Ray &v) = default;
    Ray &operator=(const Ray &v) = default;

    Vec3 operator()(double t) const{
        Vec3 res_vec = org + dir * t;
        return res_vec;
    }


    Vec3 org;
    Vec3 dir;
    Vec3::value_type tnear;
    Vec3::value_type tfar;
};


/************************************************    Implementation    ************************************************/

inline Ray::Ray()               : org(),
                                  dir(),
                                  tnear(std::numeric_limits<Vec3::value_type>::min()),
                                  tfar(std::numeric_limits<Vec3::value_type>::max()) { }

inline Ray::Ray(Vec3 o, Vec3 d) : org(o),
                                  dir(d),
                                  tnear(std::numeric_limits<Vec3::value_type>::min()),
                                  tfar(std::numeric_limits<Vec3::value_type>::max()) { }


/*
 * Functions for copying data from Ray to RTC Ray structures and vice versa.
 * RTCRayHit is used for normal ray intersections. Now contains two sub-structures ray and hit (compared to Embree2).
 * RTCRay is used for occlusion ray intersections.
 */
inline void Ray_to_RTCRayHit(const Ray& ray, RTCRayHit& rtc_ray, unsigned int geomID = (int)RTC_INVALID_GEOMETRY_ID,
                        unsigned int primID = (int)RTC_INVALID_GEOMETRY_ID, unsigned int mask = -1, float time = 0.0f) {
    rtc_ray.ray.org_x = ray.org.x;
    rtc_ray.ray.org_y = ray.org.y;
    rtc_ray.ray.org_z = ray.org.z;

    rtc_ray.ray.dir_x = ray.dir.x;
    rtc_ray.ray.dir_y = ray.dir.y;
    rtc_ray.ray.dir_z = ray.dir.z;

    rtc_ray.ray.tnear = ray.tnear;
    rtc_ray.ray.tfar = ray.tfar;

    rtc_ray.hit.geomID = geomID;
    rtc_ray.hit.primID = primID;
    rtc_ray.ray.mask = mask;
    rtc_ray.ray.time = time;
}

inline void Ray_to_RTCRay(const Ray& ray, RTCRay& rtc_ray, unsigned int mask = -1, float time = 0.0f) {
    rtc_ray.org_x = ray.org.x;
    rtc_ray.org_y = ray.org.y;
    rtc_ray.org_z = ray.org.z;

    rtc_ray.dir_x = ray.dir.x;
    rtc_ray.dir_y = ray.dir.y;
    rtc_ray.dir_z = ray.dir.z;

    rtc_ray.tnear = ray.tnear;
    rtc_ray.tfar = ray.tfar;

    rtc_ray.mask = mask;
    rtc_ray.time = time;
}

inline void RTCRayHit_to_Ray(const RTCRayHit& rtcray, Ray& ray) {
    ray.org = Vec3(rtcray.ray.org_x, rtcray.ray.org_y, rtcray.ray.org_z);
    ray.dir = Vec3(rtcray.ray.dir_x, rtcray.ray.dir_y, rtcray.ray.dir_z);
    ray.tfar = rtcray.ray.tfar;
    ray.tnear = rtcray.ray.tnear;
}

