#include <iostream>
#include <vector>
#include <cmath>
#include "spinor.h"
#include "cuba.h"
#include <random>
#include <numbers>

template <typename T>
Vector<T, 4> sphericalMap(T t[]) {
    auto constexpr pi = M_PI;
    T r = 1.0 / t[0] - 1;
    T c1 = ::cos(pi*t[1]);
    T s1 = ::sin(pi*t[1]);
    T c2 = ::cos(pi*t[2]);
    T s2 = ::sin(pi*t[2]);
    T c3 = ::cos(2*pi*t[3]);
    T s3 = ::sin(2*pi*t[3]);
    return Vector<T, 4>{r*c1, r*s1*c2, r*s1*s2*c3, r*s1*s2*s3};
}

template <typename T>
T sphericalMapScale(T t[]) {
    auto constexpr pi = M_PI;
    T r = 1.0 / t[0] - 1;
    T s1 = ::sin(pi*t[1]);
    T s2 = ::sin(pi*t[2]);
    return ::pow(r, 3) * s1 * s1 * s2 * 2 * ::pow(pi, 3) / t[0] / t[0];
}

template <typename T>
Vector<T, 4> sphericalMap3D(T t[]) {
    auto constexpr pi = M_PI;
    T r = 1.0 / t[1] - 1;
    T c1 = ::cos(pi*t[2]);
    T s1 = ::sin(pi*t[2]);
    T c2 = ::cos(2*pi*t[3]);
    T s2 = ::sin(2*pi*t[3]);
    return Vector<T, 4>{t[0], r*c1, r*s1*c2, r*s1*s2};
}

template <typename T>
T sphericalMapScale3D(T t[]) {
    auto constexpr pi = M_PI;
    T r = 1.0 / t[1] - 1;
    T s1 = ::sin(pi*t[2]);
    return r * r * s1 * 2 * ::pow(pi, 2) / t[1] / t[1];
}

template <typename T>
auto dot(T v1, T v2) {
    return v1[0] * v2[0] - v1[1] * v2[1] - v1[2] * v2[2] - v1[3] * v2[3];
}

template <typename T>
auto square(T v) {
    return dot(v, v);
}

template <typename T>
auto metric(T v) {
    return T{v[0], -v[1], -v[2], -v[3]};
}

template <typename T, typename U>
auto convert(U u) {
    return T{u[0], u[1], u[2], u[3]};
}

std::array<Vector<std::complex<double>, 4>, 4> gammaUp = {
    Vector<std::complex<double>, 4>{1.0, 0, 0, 0},
    Vector<std::complex<double>, 4>{0.0, 1, 0, 0},
    Vector<std::complex<double>, 4>{0.0, 0, 1, 0},
    Vector<std::complex<double>, 4>{0.0, 0, 0, 1},
};

std::array<Vector<std::complex<double>, 4>, 4> gammaDown = {
    Vector<std::complex<double>, 4>{1.0, 0, 0, 0},
    Vector<std::complex<double>, 4>{0.0, -1, 0, 0},
    Vector<std::complex<double>, 4>{0.0, 0, -1, 0},
    Vector<std::complex<double>, 4>{0.0, 0, 0, -1},
};
