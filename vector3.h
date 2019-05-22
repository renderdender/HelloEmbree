#pragma once

#include <iostream>
#include <cstddef>
#include <cmath>
#include <cassert>


/*************************************************    Declaration    **************************************************/

template <typename T>
class Vector3T {

public:
    static const std::size_t dim = 3;

    // Initialize with default value
    explicit Vector3T();
    // Initialize with a certain value
    explicit Vector3T(T val);
    // Initialize with 3 values
    explicit Vector3T(T val_x, T val_y, T val_z);
    // Initialize with an array of the same type
    explicit Vector3T(T val_array[3]);

    ~Vector3T() = default;

    // Copy constructor with Vector3T of type K
    template <typename K>
    Vector3T(const Vector3T<K> &v);
    // Default copy constructor and assignment operator
    Vector3T(const Vector3T<T> &v) = default;
    Vector3T &operator=(const Vector3T &v) = default;


    T& operator[](const size_t i);
    T operator[](const size_t i) const;

    typedef T value_type;
    // Vector components
    T x, y, z;
};

// Basic vector arithmetics
template <typename T> Vector3T<T> operator+ (const Vector3T<T>& lhs, const Vector3T<T>& rhs);
template <typename T> Vector3T<T> operator+=(      Vector3T<T>& lhs, const Vector3T<T>& rhs);

template <typename T> Vector3T<T> operator- (const Vector3T<T>& lhs);
template <typename T> Vector3T<T> operator- (const Vector3T<T>& lhs, const Vector3T<T>& rhs);
template <typename T> Vector3T<T> operator-=(      Vector3T<T>& lhs, const Vector3T<T>& rhs);

template <typename T> Vector3T<T> operator* (const Vector3T<T>& lhs, const T rhs);
template <typename T> Vector3T<T> operator* (const T lhs,            const Vector3T<T>& rhs);
template <typename T> Vector3T<T> operator* (const Vector3T<T>& lhs, const Vector3T<T>& rhs);
template <typename T> Vector3T<T> operator*=(      Vector3T<T>& lhs, const T rhs);
template <typename T> Vector3T<T> operator*=(      Vector3T<T>& lhs, const Vector3T<T>& rhs);

template <typename T> Vector3T<T> operator/ (const Vector3T<T>& lhs, const T rhs);
template <typename T> Vector3T<T> operator/ (const Vector3T<T>& lhs, const Vector3T<T>& rhs);
template <typename T> Vector3T<T> operator/=(      Vector3T<T>& lhs, const T rhs);
template <typename T> Vector3T<T> operator/=(      Vector3T<T>& lhs, const Vector3T<T>& rhs);

// Equality / Inequality tests
template <typename T> bool        operator!=(const Vector3T<T>& lhs, const Vector3T<T>& rhs);
template <typename T> bool        operator==(const Vector3T<T>& lhs, const Vector3T<T>& rhs);

// Dot product, Cross product, Basis change
template <typename T> T           dot       (const Vector3T<T>& lhs, const Vector3T<T>& rhs);
template <typename T> Vector3T<T> cross     (const Vector3T<T>& lhs, const Vector3T<T>& rhs);
template <typename T> Vector3T<T> basis     (const Vector3T<T>& lhs, const Vector3T<T>& rhs);

// Norms
template <typename T> T           norm2     (const Vector3T<T>& lhs);
template <typename T> T           norm      (const Vector3T<T>& lhs);
template <typename T> Vector3T<T> normalize (      Vector3T<T>  lhs);

/********************************************  Typedefs for common vectors  *******************************************/

typedef Vector3T<unsigned int>  Vector3I;
typedef Vector3T<float>         Vector3F;
typedef Vector3T<double>        Vector3D;
typedef Vector3I                Vec3i;
typedef Vector3F                Vec3f;
typedef Vector3D                Vec3;

/***************************************************  Implementation  *************************************************/
template <typename T>
inline Vector3T<T>::Vector3T()                            : x(std::numeric_limits<T>::min()),
                                                            y(std::numeric_limits<T>::min()),
                                                            z(std::numeric_limits<T>::min()) { }

template <typename T>
inline Vector3T<T>::Vector3T(T val)                       : x(val),
                                                            y(val),
                                                            z(val) { }
template <typename T>
inline Vector3T<T>::Vector3T(T val_x, T val_y, T val_z)   : x(val_x),
                                                            y(val_y),
                                                            z(val_z) { }

template <typename T>
inline Vector3T<T>::Vector3T(T val_array[3])              : x(val_array[0]),
                                                            y(val_array[1]),
                                                            z(val_array[2]) { }

template <typename T>
template <typename K>
inline Vector3T<T>::Vector3T(const Vector3T<K> &v)        : x(static_cast<T>(v.x)),
                                                            y(static_cast<T>(v.y)),
                                                            z(static_cast<T>(v.z)) { }

template <typename T>
inline T& Vector3T<T>::operator[](const size_t i){
    assert(i < dim);

    return (&x)[i];
}

template <typename T>
inline T Vector3T<T>::operator[](const size_t i) const{
    assert(i < dim);

    return (&x)[i];
}

// Basic vector arithmetics
template <typename T>
inline Vector3T<T> operator+ (const Vector3T<T>& lhs, const Vector3T<T>& rhs){
    Vector3T<T> res_vec(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);

    return res_vec;
}

template <typename T>
inline Vector3T<T> operator+=(      Vector3T<T>& lhs, const Vector3T<T>& rhs){
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;

    return lhs;
}

template <typename T>
inline Vector3T<T> operator- (const Vector3T<T>& lhs){
    Vector3T<T> res_vec(-lhs.x, -lhs.y, -lhs.z);

    return res_vec;
}

template <typename T>
inline Vector3T<T> operator- (const Vector3T<T>& lhs, const Vector3T<T>& rhs){
    Vector3T<T> res_vec(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);

    return res_vec;
}

template <typename T>
inline Vector3T<T> operator-=(      Vector3T<T>& lhs, const Vector3T<T>& rhs){
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;

    return lhs;
}

template <typename T>
inline Vector3T<T> operator* (const Vector3T<T>& lhs, const T rhs){
    Vector3T<T> res_vec(lhs.x * rhs, lhs.y * rhs, lhs.z * rhs);

    return res_vec;
}

template <typename T>
inline Vector3T<T> operator* (const T lhs,           const Vector3T<T>& rhs){
    Vector3T<T> res_vec(rhs.x * lhs, rhs.y * lhs, rhs.z * lhs);

    return res_vec;
}

template <typename T>
inline Vector3T<T> operator* (const Vector3T<T>& lhs, const Vector3T<T>& rhs){
    Vector3T<T> res_vec(lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z);

    return res_vec;
}

template <typename T>
inline Vector3T<T> operator*=(      Vector3T<T>& lhs, const T rhs){
    lhs.x *= rhs;
    lhs.y *= rhs;
    lhs.z *= rhs;

    return lhs;
}

template <typename T>
inline Vector3T<T> operator*=(      Vector3T<T>& lhs, const Vector3T<T>& rhs){
    lhs.x *= rhs.x;
    lhs.y *= rhs.y;
    lhs.z *= rhs.z;

    return lhs;
}

template <typename T>
inline Vector3T<T> operator/ (const Vector3T<T>& lhs, const T rhs){
    Vector3T<T> res_vec(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs);

    return res_vec;
}

template <typename T>
inline Vector3T<T> operator/ (const Vector3T<T>& lhs, const Vector3T<T>& rhs){
    Vector3T<T> res_vec(lhs.x / rhs.x, lhs.y / rhs.y, lhs.z / rhs.z);

    return res_vec;
}

template <typename T>
inline Vector3T<T> operator/=(      Vector3T<T>& lhs, const T rhs){
    lhs.x /= rhs;
    lhs.y /= rhs;
    lhs.z /= rhs;

    return lhs;
}

template <typename T>
inline Vector3T<T> operator/=(      Vector3T<T>& lhs, const Vector3T<T>& rhs){
    lhs.x /= rhs.x;
    lhs.y /= rhs.y;
    lhs.z /= rhs.z;

    return lhs;
}

// Equality / Inequality tests
template <typename T>
inline bool        operator!=(const Vector3T<T>& lhs, const Vector3T<T>& rhs){

    return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z);
}

template <typename T>
inline bool        operator==(const Vector3T<T>& lhs, const Vector3T<T>& rhs){

    return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

// Dot product, Cross product.
template <typename T>
inline T           dot(const Vector3T<T>& lhs, const Vector3T<T>& rhs){
    T res_val(lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z);

    return res_val;
}

template <typename T>
inline Vector3T<T> cross(const Vector3T<T>& lhs, const Vector3T<T>& rhs){
    Vector3T<T> res_vec(lhs.y * rhs.z - rhs.y * lhs.z,
                       lhs.z * rhs.x - rhs.z * lhs.x,
                       lhs.x * rhs.y - rhs.x * lhs.y);

    return res_vec;
}

template <typename T>
inline Vector3T<T> basis(const Vector3T<T> &lhs, const Vector3T<T> &rhs){
    Vec3 u, w, v = rhs;
    if (rhs.z < -0.9999999) {
        u = Vec3(0.0, -1.0, 0.0);
        w = Vec3(-1.0, 0.0, 0.0);
    }
    else {
        const double a = 1.0 / (1.0 + rhs.z);
        const double b = -rhs.x * rhs.y * a;
        u = Vec3(1.0 - rhs.x * rhs.x * a, b, -rhs.x);
        w = Vec3(b, 1.0 - rhs.y * rhs.y * a, -rhs.y);
    }
    return Vec3(dot(lhs, Vec3(u.x, v.x, w.x)), dot(lhs, Vec3(u.y, v.y, w.y)), dot(lhs, Vec3(u.z, v.z, w.z)));
}

// Norms
template <typename T>
inline T           norm2(const Vector3T<T>& lhs){

    return dot(lhs, lhs);
}

template <typename T>
inline T           norm(const Vector3T<T>& lhs){

    return std::sqrt(norm2(lhs));
}

template <typename T>
inline Vector3T<T> normalize(Vector3T<T> lhs){
    T div = 1.0 / norm(lhs);
    lhs.x *= div;
    lhs.y *= div;
    lhs.z *= div;

    return lhs;
}

