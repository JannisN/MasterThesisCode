#pragma once

#include "matrix.h"
#include <complex>
#include <cmath>
#include <iostream>

template <typename T, typename U>
auto createPauli(U u) {
    return Matrix<std::complex<T>, 2, 2>(u[0] + u[3], u[1] - std::complex<T>(1i) * u[2], u[1] + std::complex<T>(1i) * u[2], u[0] - u[3]);
}

template <typename T, typename U>
auto createAntiPauli(U u) {
    return Matrix<std::complex<T>, 2, 2>(u[0] - u[3], -u[1] + std::complex<T>(1i) * u[2], -u[1] - std::complex<T>(1i) * u[2], u[0] + u[3]);
}

template <typename T, typename U>
auto createSpinorMatrix(U u) {
    T a0 = std::sqrt((u[0] + std::sqrt(u[0] * u[0] - u[1] * u[1] - u[2] * u[2] - u[3] * u[3])) * 0.5);
    // - signs because we have to pull the index down
    T a1 = -u[1] / (2 * a0);
    T a2 = -u[2] / (2 * a0);
    T a3 = -u[3] / (2 * a0);
    return createPauli<T>(std::array<T, 4>{a0, a1, a2, a3});
}

template <typename T, typename U>
auto createAntiSpinorMatrix(U u) {
    T a0 = std::sqrt((u[0] + std::sqrt(u[0] * u[0] - u[1] * u[1] - u[2] * u[2] - u[3] * u[3])) * 0.5);
    T a1 = -u[1] / (2 * a0);
    T a2 = -u[2] / (2 * a0);
    T a3 = -u[3] / (2 * a0);
    return createAntiPauli<T>(std::array<T, 4>{a0, a1, a2, a3});
}

template <typename T, typename A>
auto createGammaChainFromVectorsMassless(A vectors, uint vectorCount) {
    Matrix<std::complex<T>, 2, 2> m1 = Matrix<std::complex<T>, 2, 2>::Identity();
    Matrix<std::complex<T>, 2, 2> m2 = Matrix<std::complex<T>, 2, 2>::Identity();
    for (uint i = 0; i < vectorCount; i++) {
        // if i is even
        if (i % 2 == 0) {
            m1 *= createPauli<T>(vectors[i]);
            m2 *= createAntiPauli<T>(vectors[i]);
        } else {
            m2 *= createPauli<T>(vectors[i]);
            m1 *= createAntiPauli<T>(vectors[i]);
        }
    }
    if (vectorCount % 2 == 0) {
        return Matrix<std::complex<T>, 4, 4> (
            m1[0], m1[1], 0, 0,
            m1[2], m1[3], 0, 0,
            0, 0, m2[0], m2[1],
            0, 0, m2[2], m2[3]
        );
    } else {
        return Matrix<std::complex<T>, 4, 4> (
            0, 0, m1[0], m1[1],
            0, 0, m1[2], m1[3],
            m2[0], m2[1], 0, 0,
            m2[2], m2[3], 0, 0
        );
    }
}

template <typename T, typename A, typename B>
auto createGammaChainFromVectorsMassive(A vectors, uint vectorCount, B masses) {
    Matrix<Matrix<std::complex<T>, 2, 2>, 2, 2> m(Matrix<std::complex<T>, 2, 2>::Identity(), Matrix<std::complex<T>, 2, 2>(), Matrix<std::complex<T>, 2, 2>(), Matrix<std::complex<T>, 2, 2>::Identity());
    for (uint i = 0; i < vectorCount; i++) {
        m *= Matrix<Matrix<std::complex<T>, 2, 2>, 2, 2>(
            Matrix<std::complex<T>, 2, 2>(masses[i], 0, 0, masses[i]), createPauli<T>(vectors[i]),
            createAntiPauli<T>(vectors[i]), Matrix<std::complex<T>, 2, 2>(masses[i], 0, 0, masses[i])
        );
    }
    return Matrix<std::complex<T>, 4, 4> (
        m[0][0], m[0][1], m[1][0], m[1][1],
        m[0][2], m[0][3], m[1][2], m[1][3],
        m[2][0], m[2][1], m[3][0], m[3][1],
        m[2][2], m[2][3], m[3][2], m[3][3]
    );
}

template <typename T, typename A, typename Momentum, typename Spin>
auto computeSpinorProductMassless(bool isVbarU, Momentum momentumLeft, Momentum momentumRight, Spin spinLeft, Spin spinRight, A vectors, uint vectorCount) {
    auto m = createGammaChainFromVectorsMassless<T>(vectors, vectorCount);
    if (isVbarU) {
        auto u1 = createSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto u2 = createAntiSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto v1 = createSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        auto v2 = -createAntiSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        // debug
        //std::cout << u1[0] << " " << u1[1] << " " << u2[0] << " " << u2[1] << std::endl;
        //std::cout << v1[0] << " " << v1[1] << " " << v2[0] << " " << v2[1] << std::endl;
        Vector<std::complex<T>, 4> u(u1[0], u1[1], u2[0], u2[1]);
        // first v2 then v1 due to the 0th gamma matrix
        Vector<std::complex<T>, 4> v(v2[0], v2[1], v1[0], v1[1]);
        return (adjoint(v) * m * u)[0];
    } else {
        auto v1 = createSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto v2 = -createAntiSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto u1 = createSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        auto u2 = createAntiSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        Vector<std::complex<T>, 4> u(u2[0], u2[1], u1[0], u1[1]);
        Vector<std::complex<T>, 4> v(v1[0], v1[1], v2[0], v2[1]);
        return (adjoint(u) * m * v)[0];
    }
}

template <typename T, typename A, typename B, typename Momentum, typename Spin>
auto computeSpinorProductMassive(bool isVbarU, Momentum momentumLeft, Momentum momentumRight, Spin spinLeft, Spin spinRight, A vectors, uint vectorCount, B masses) {
    auto m = createGammaChainFromVectorsMassive<T>(vectors, vectorCount, masses);
    if (isVbarU) {
        auto u1 = createSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto u2 = createAntiSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto v1 = createSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        auto v2 = -createAntiSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        Vector<std::complex<T>, 4> u(u1[0], u1[1], u2[0], u2[1]);
        // first v2 then v1 due to the 0th gamma matrix
        Vector<std::complex<T>, 4> v(v2[0], v2[1], v1[0], v1[1]);
        return (adjoint(v) * m * u)[0];
    } else {
        auto v1 = createSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto v2 = -createAntiSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto u1 = createSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        auto u2 = createAntiSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        Vector<std::complex<T>, 4> u(u2[0], u2[1], u1[0], u1[1]);
        Vector<std::complex<T>, 4> v(v1[0], v1[1], v2[0], v2[1]);
        return (adjoint(u) * m * v)[0];
    }
}

template <typename T, typename A>
auto createGammaChainFromVectorsMasslessChiral(A vectors, uint vectorCount) {
    Matrix<std::complex<T>, 2, 2> m1 = Matrix<std::complex<T>, 2, 2>::Identity();
    Matrix<std::complex<T>, 2, 2> m2 = Matrix<std::complex<T>, 2, 2>::Identity();
    // because the 5th gamma matrix is diagonal
    uint chiralCount = 0;
    for (uint i = 0; i < vectorCount; i++) {
        if ((i + chiralCount) % 2 == 0) {
            m1 *= createPauli<T>(vectors[i]) + Matrix<std::complex<T>, 2, 2>(-vectors[i][4], 0, 0, -vectors[i][4]);
            m2 *= createAntiPauli<T>(vectors[i]) + Matrix<std::complex<T>, 2, 2>(vectors[i][4], 0, 0, vectors[i][4]);
        } else {
            m2 *= createPauli<T>(vectors[i]) + Matrix<std::complex<T>, 2, 2>(-vectors[i][4], 0, 0, -vectors[i][4]);
            m1 *= createAntiPauli<T>(vectors[i]) + Matrix<std::complex<T>, 2, 2>(vectors[i][4], 0, 0, vectors[i][4]);
        }
        if (vectors[i][4] != 0.0) {
            chiralCount++;
        }
    }
    if ((vectorCount + chiralCount) % 2 == 0) {
        return Matrix<std::complex<T>, 4, 4> (
            m1[0], m1[1], 0, 0,
            m1[2], m1[3], 0, 0,
            0, 0, m2[0], m2[1],
            0, 0, m2[2], m2[3]
        );
    } else {
        return Matrix<std::complex<T>, 4, 4> (
            0, 0, m1[0], m1[1],
            0, 0, m1[2], m1[3],
            m2[0], m2[1], 0, 0,
            m2[2], m2[3], 0, 0
        );
    }
}

template <typename T, typename A, typename B>
auto createGammaChainFromVectorsMassiveChiral(A vectors, uint vectorCount, B masses) {
    Matrix<Matrix<std::complex<T>, 2, 2>, 2, 2> m(Matrix<std::complex<T>, 2, 2>::Identity(), Matrix<std::complex<T>, 2, 2>(), Matrix<std::complex<T>, 2, 2>(), Matrix<std::complex<T>, 2, 2>::Identity());
    for (uint i = 0; i < vectorCount; i++) {
        m *= Matrix<Matrix<std::complex<T>, 2, 2>, 2, 2>(
            Matrix<std::complex<T>, 2, 2>(masses[i] - vectors[i][4], 0, 0, masses[i] - vectors[i][4]), createPauli<T>(vectors[i]),
            createAntiPauli<T>(vectors[i]), Matrix<std::complex<T>, 2, 2>(masses[i] + vectors[i][4], 0, 0, masses[i] + vectors[i][4])
        );
    }
    return Matrix<std::complex<T>, 4, 4> (
        m[0][0], m[0][1], m[1][0], m[1][1],
        m[0][2], m[0][3], m[1][2], m[1][3],
        m[2][0], m[2][1], m[3][0], m[3][1],
        m[2][2], m[2][3], m[3][2], m[3][3]
    );
}

template <typename T, typename A, typename Momentum, typename Spin>
auto computeSpinorProductMasslessChiral(bool isVbarU, Momentum momentumLeft, Momentum momentumRight, Spin spinLeft, Spin spinRight, A vectors, uint vectorCount) {
    auto m = createGammaChainFromVectorsMasslessChiral<T>(vectors, vectorCount);
    if (isVbarU) {
        auto u1 = createSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto u2 = createAntiSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto v1 = createSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        auto v2 = -createAntiSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        // debug
        //std::cout << u1[0] << " " << u1[1] << " " << u2[0] << " " << u2[1] << std::endl;
        //std::cout << v1[0] << " " << v1[1] << " " << v2[0] << " " << v2[1] << std::endl;
        Vector<std::complex<T>, 4> u(u1[0], u1[1], u2[0], u2[1]);
        // first v2 then v1 due to the 0th gamma matrix
        Vector<std::complex<T>, 4> v(v2[0], v2[1], v1[0], v1[1]);
        return (adjoint(v) * m * u)[0];
    } else {
        auto v1 = createSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto v2 = -createAntiSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto u1 = createSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        auto u2 = createAntiSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        Vector<std::complex<T>, 4> u(u2[0], u2[1], u1[0], u1[1]);
        Vector<std::complex<T>, 4> v(v1[0], v1[1], v2[0], v2[1]);
        return (adjoint(u) * m * v)[0];
    }
}

template <typename T, typename A, typename B, typename Momentum, typename Spin>
auto computeSpinorProductMassiveChiral(bool isVbarU, Momentum momentumLeft, Momentum momentumRight, Spin spinLeft, Spin spinRight, A vectors, uint vectorCount, B masses) {
    auto m = createGammaChainFromVectorsMassiveChiral<T>(vectors, vectorCount, masses);
    if (isVbarU) {
        auto u1 = createSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto u2 = createAntiSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto v1 = createSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        auto v2 = -createAntiSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        Vector<std::complex<T>, 4> u(u1[0], u1[1], u2[0], u2[1]);
        // first v2 then v1 due to the 0th gamma matrix
        Vector<std::complex<T>, 4> v(v2[0], v2[1], v1[0], v1[1]);
        return (adjoint(v) * m * u)[0];
    } else {
        auto v1 = createSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto v2 = -createAntiSpinorMatrix<T>(momentumRight) * Vector<std::complex<T>, 2>(spinRight[0], spinRight[1]);
        auto u1 = createSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        auto u2 = createAntiSpinorMatrix<T>(momentumLeft) * Vector<std::complex<T>, 2>(spinLeft[0], spinLeft[1]);
        Vector<std::complex<T>, 4> u(u2[0], u2[1], u1[0], u1[1]);
        Vector<std::complex<T>, 4> v(v1[0], v1[1], v2[0], v2[1]);
        return (adjoint(u) * m * v)[0];
    }
}
