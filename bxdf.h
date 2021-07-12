#pragma once

#include "vector3.h"
#include "utils.h"

enum BxDF_TYPE {
    DIFFUSE,
    SPECULAR,
    REFRACTIVE
};


// Ideal specular reflection and it's discreet probability
// 'direction' is oriented TO the surface
// in other renderers (e.g. Mitsuba/pbrt) it is the opposite (= '2*n*dot(d,n) - d')
inline Vec3 reflect_at_normal(const Vec3 &dir, const Vec3 &n){
    return dir - n * 2.0 * dot(n, dir);
}

// Just a 'typedef' for convenience
// 'direction' is oriented AWAY FROM the surface
inline Vec3 specular_reflection(const Vec3 &wi, const Vec3 &n){
    return reflect_at_normal(-wi, n);
}

// Uniform sample on the sphere and it's PDF
inline Vec3 uniform_sphere(const double rnd1, const double rnd2){
    const double cos_theta = 1.0 - 2.0 * rnd1;
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    const double phi = 2.0 * D_PI * rnd2;
    Vec3 res_vec = Vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);

    return res_vec;
}

inline double uniform_sphere_pdf() {
    return D_INV4_PI;
}

// Uniform sample on the hemisphere and it's PDF
inline Vec3 uniform_hemisphere(const double rnd1, const double rnd2){
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - rnd1 * rnd1));
    const double phi = 2.0 * D_PI * rnd2;
    Vec3 res_vec = Vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, rnd2);

    return res_vec;
}

inline double uniform_hemisphere_pdf() {
    return D_INV2_PI;
}

// Uniform sample on the cosine weighted hemisphere and it's PDF
inline Vec3 cosine_weighted_hemisphere(const double rnd1, const double rnd2){
    const double cos_theta = std::sqrt(1.0 - rnd1);
    const double sin_theta = std::sqrt(rnd1);
    const double phi = 2.0 * D_PI * rnd2;
    return Vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
}

inline double cosine_wighted_hemisphere_pdf(double cos_theta) {
    return cos_theta * D_INV_PI;
}

inline double cosine_wighted_hemisphere_pdf(const Vec3 wi, const Vec3 n, const Vec3 wo){
    double cos_theta = std::abs(dot(wo, n));
    return cos_theta * D_INV_PI;
}

