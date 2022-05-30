#include "integration.h"

template <typename T>
struct ProcessInfo {
    Vector<T, 4> p1, p2, q1, q2;
    Vector<std::complex<T>, 2> spinP1, spinP2, spinQ1, spinQ2;
    T mass;
    T uvScale;
    T t1, t2;
    T ca, cf;
};

template <std::complex<double> F(Vector<std::complex<double>, 4>, ProcessInfo<double>, double, std::complex<double>)>
int integrateDiagramDouble(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale = sphericalMapScale(x);
    Vector<double, 4> kReal = sphericalMap<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]), kReal[1], kReal[2], kReal[3]};
    double delta = 0.0;
    std::complex<double> direction = std::complex<double>(1.0i);
    std::complex<double> result = std::complex<double>(1.0i) * F(k, *static_cast<ProcessInfo<double>*>(userdata), delta, direction) * scale;
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

template <typename T>
std::complex<T> masslessTriangleDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    //std::cout << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> kReal = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p2) + k) + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = 0.0;
    std::complex<T> uvNumerator = 0.0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 5> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(convert<Vector<std::complex<T>, 4>>(boxInfo.p2) + kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[j]), metric(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            std::array<Vector<std::complex<T>, 4>, 1> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[j]) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 5)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
            std::array<Vector<std::complex<T>, 4>, 5> momentaMasslessUV = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[j]), metric(- kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            uvNumerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMasslessUV, 5)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
        }
    }
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> counterIR = 2.0 * (
        2 * dot(boxInfo.p1, boxInfo.p2) / (((square(k) + std::complex<T>(1.0i) * Delta) * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta) * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p2) + k))) + std::complex<T>(1.0i) * Delta)
        - 1.0 / (square(k) * square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta)
        - 1.0 / (square(k) * square(convert<Vector<std::complex<T>, 4>>(boxInfo.p2) + k) + std::complex<T>(1.0i) * Delta)
        + 2.0 / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta)
    );
    std::complex<T> counterUV1 = 2.0 * square(kReal) / (square(k) * boxInfo.uvScale * boxInfo.uvScale) / (square(k) * boxInfo.uvScale * boxInfo.uvScale) / (square(k) * boxInfo.uvScale * boxInfo.uvScale);
    std::complex<T> counterUV2 = -4.0 * computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, std::array<Vector<std::complex<T>, 4>, 1>{metric(kReal)}, 1)
        * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, std::array<Vector<std::complex<T>, 4>, 1>{metric(kReal)}, 1)
        / (square(k) * boxInfo.uvScale * boxInfo.uvScale) / (square(k) * boxInfo.uvScale * boxInfo.uvScale) / (square(k) * boxInfo.uvScale * boxInfo.uvScale);

    std::complex<T> tree = 0;
    for (int i = 0; i < 4; i++) {
        std::array<Vector<std::complex<T>, 4>, 1> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
        std::array<Vector<std::complex<T>, 4>, 1> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
        tree += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
            * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
    }
    counterIR *= tree;
    counterUV1 *= tree;
    return (numerator / denominator) - counterIR - uvNumerator / uvDenominator;
    //return (numerator / denominator) - counterIR - counterUV1 - counterUV2;
}

template <typename T>
std::complex<T> massiveTriangleDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> kReal = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    //std::cout << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << std::endl;
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q2) - k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q1) + k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = 0.0;
    std::complex<T> uvNumerator = 0.0;
    std::array<T, 5> masses = {0, boxInfo.mass, 0, boxInfo.mass, 0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 5> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(convert<Vector<std::complex<T>, 4>>(boxInfo.q1) + kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[j]), metric(-convert<Vector<std::complex<T>, 4>>(boxInfo.q2) + kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            std::array<Vector<std::complex<T>, 4>, 1> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaDown[j]) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
                * computeSpinorProductMassive<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 5, masses);
            std::array<Vector<std::complex<T>, 4>, 5> momentaMassiveUV = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[j]), metric(kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            uvNumerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassiveUV, 5);
        }
    }
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta)
        * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta)
        * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> counterIR = 2.0 * (
        -2 * dot(boxInfo.q1, boxInfo.q2) / (((square(k) + std::complex<T>(1.0i) * Delta) * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q2) - k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta) * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q1) + k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta)))
    );
    std::complex<T> counterUV1 = 2.0 * square(kReal) / (square(k) * boxInfo.uvScale * boxInfo.uvScale) / (square(k) * boxInfo.uvScale * boxInfo.uvScale) / (square(k) * boxInfo.uvScale * boxInfo.uvScale);
    std::complex<T> counterUV2 = -4.0 * computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, std::array<Vector<std::complex<T>, 4>, 1>{metric(kReal)}, 1)
        * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, std::array<Vector<std::complex<T>, 4>, 1>{metric(kReal)}, 1)
        / (square(k) * boxInfo.uvScale * boxInfo.uvScale) / (square(k) * boxInfo.uvScale * boxInfo.uvScale) / (square(k) * boxInfo.uvScale * boxInfo.uvScale);

    std::complex<T> tree = 0;
    for (int i = 0; i < 4; i++) {
        std::array<Vector<std::complex<T>, 4>, 1> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
        std::array<Vector<std::complex<T>, 4>, 1> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
        tree += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
            * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
    }
    counterIR *= tree;
    counterUV1 *= tree;
    //std::cout << (denominator) <<  "     " << uvDenominator << std::endl;
    return (numerator / denominator) - counterIR - uvNumerator / uvDenominator;
    //return (numerator / denominator) - counterIR - counterUV1 - counterUV2;
}

template <typename T>
std::complex<T> masslessGluonTriangleDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> kReal = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p2) + k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta);
    Vector<std::complex<T>, 4> vRho = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - convert<Vector<std::complex<T>, 4>>(boxInfo.p2) - kReal - kReal;
    Vector<std::complex<T>, 4> vNu = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2) + kReal;
    Vector<std::complex<T>, 4> vMu = -convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - convert<Vector<std::complex<T>, 4>>(boxInfo.p2) + kReal;
    Vector<std::complex<T>, 4> vRhoUV = - kReal - kReal;
    Vector<std::complex<T>, 4> vNuUV = kReal;
    Vector<std::complex<T>, 4> vMuUV = kReal;
    std::complex<T> numerator = 0.0;
    std::complex<T> uvNumerator = 0.0;
    for (int i = 0; i < 4; i++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 1> momentaMassive = { metric(vRho) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 3)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
            std::array<Vector<std::complex<T>, 4>, 1> momentaMassiveUV = { metric(vRhoUV) };
            uvNumerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 3)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassiveUV, 1);
    }
    for (int i = 0; i < 4; i++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 1> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(kReal), metric(vNu) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMasslessUV = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(kReal), metric(vNuUV) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 3)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
            uvNumerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMasslessUV, 3)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
    }
    for (int i = 0; i < 4; i++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 1> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassless = { metric(vNu), metric(kReal), convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMasslessUV = { metric(vNuUV), metric(kReal), convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 3)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
            uvNumerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMasslessUV, 3)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
    }

    std::complex<T> counter = uvNumerator / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> tree = 0;
    for (int i = 0; i < 4; i++) {
        std::array<Vector<std::complex<T>, 4>, 3> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
        std::array<Vector<std::complex<T>, 4>, 3> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
        tree += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
            * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
    }
    counter += tree * (-2.0 / ((square(k) + std::complex<T>(1.0i) * Delta) * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta))
                       -2.0 / ((square(k) + std::complex<T>(1.0i) * Delta) * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p2) + k) + std::complex<T>(1.0i) * Delta))
                       +4.0 / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta)
                       );
    return (numerator / denominator) - counter;
}

template <typename T>
std::complex<T> massiveGluonTriangleDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> kReal = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q1) + k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q2) - k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta);
    Vector<std::complex<T>, 4> vRho = convert<Vector<std::complex<T>, 4>>(boxInfo.q2) - convert<Vector<std::complex<T>, 4>>(boxInfo.q1) - kReal - kReal;
    Vector<std::complex<T>, 4> vNu = convert<Vector<std::complex<T>, 4>>(boxInfo.q2) + convert<Vector<std::complex<T>, 4>>(boxInfo.q1) + convert<Vector<std::complex<T>, 4>>(boxInfo.q1) + kReal;
    Vector<std::complex<T>, 4> vMu = -convert<Vector<std::complex<T>, 4>>(boxInfo.q2) - convert<Vector<std::complex<T>, 4>>(boxInfo.q2) - convert<Vector<std::complex<T>, 4>>(boxInfo.q1) + kReal;
    Vector<std::complex<T>, 4> vRhoUV = - kReal - kReal;
    Vector<std::complex<T>, 4> vNuUV = kReal;
    Vector<std::complex<T>, 4> vMuUV = kReal;
    std::complex<T> numerator = 0.0;
    std::complex<T> uvNumerator = 0.0;
    std::array<T, 3> masses = {0, boxInfo.mass, 0};
    for (int i = 0; i < 4; i++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 1> momentaMassless = { metric(vRho) };
            std::array<Vector<std::complex<T>, 4>, 1> momentaMasslessUV = { metric(vRhoUV) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
                * computeSpinorProductMassive<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 3, masses);
            uvNumerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMasslessUV, 1)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 3);
    }
    for (int i = 0; i < 4; i++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 1> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(kReal), metric(vNu) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassiveUV = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(kReal), metric(vNuUV) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
                * computeSpinorProductMassive<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 3, masses);
            uvNumerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassiveUV, 3);
    }
    for (int i = 0; i < 4; i++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 1> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassive = { metric(vNu), metric(kReal), convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassiveUV = { metric(vNuUV), metric(kReal), convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
                * computeSpinorProductMassive<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 3, masses);
            uvNumerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
                * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassiveUV, 3);
    }

    std::complex<T> counter = uvNumerator / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) / (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    return (numerator / denominator) - counter;
}

template <typename T>
std::complex<T> boxDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> kReal = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1 + boxInfo.p2) - k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q1) - k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = 0.0;
    std::array<T, 3> masses = {0, boxInfo.mass, 0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]), metric(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[j]) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[j]), metric(convert<Vector<std::complex<T>, 4>>(boxInfo.q1) - kReal), convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 3)
                * computeSpinorProductMassive<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 3, masses);
        }
    }

    std::complex<T> counter = 4.0 * dot(boxInfo.p1, boxInfo.q1) / (2 * dot(boxInfo.p1, boxInfo.p2)) / (square(k) + std::complex<T>(1.0i) * Delta) / (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta)
        / (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q1) - k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta);
    counter += 4.0 * dot(boxInfo.p2, boxInfo.q2) / (2 * dot(boxInfo.p1, boxInfo.p2)) / (square(k + boxInfo.p1 + boxInfo.p2) + std::complex<T>(1.0i) * Delta) / (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta)
        / (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q1) - k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta);

    std::complex<T> tree = 0;
    for (int i = 0; i < 4; i++) {
        std::array<Vector<std::complex<T>, 4>, 3> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
        std::array<Vector<std::complex<T>, 4>, 3> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
        tree += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
            * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
    }
    counter *= tree;
    return (numerator / denominator) - counter;
}

template <typename T>
std::complex<T> crossBoxDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> kReal = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1 + boxInfo.p2) - k) + std::complex<T>(1.0i) * Delta)
        * (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q2) - k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = 0.0;
    std::array<T, 3> masses = {0, boxInfo.mass, 0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            // metric needed because momenta are assumed to have their index on the bottom with gamma matrices having their indices at the top
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]), metric(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - kReal), convert<Vector<std::complex<T>, 4>>(gammaUp[j]) };
            std::array<Vector<std::complex<T>, 4>, 3> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]), metric(kReal - convert<Vector<std::complex<T>, 4>>(boxInfo.q2)), convert<Vector<std::complex<T>, 4>>(gammaDown[j]) };
            numerator += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 3)
                * computeSpinorProductMassive<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 3, masses);
        }
    }

    std::complex<T> counter = -4.0 * dot(boxInfo.p1, boxInfo.q2) / (2 * dot(boxInfo.p1, boxInfo.p2)) / (square(k) + std::complex<T>(1.0i) * Delta) / (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta)
        / (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q2) - k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta);
    counter -= 4.0 * dot(boxInfo.p2, boxInfo.q1) / (2 * dot(boxInfo.p1, boxInfo.p2)) / (square(k + boxInfo.p1 + boxInfo.p2) + std::complex<T>(1.0i) * Delta) / (square(convert<Vector<std::complex<T>, 4>>(boxInfo.p1) - k) + std::complex<T>(1.0i) * Delta)
        / (square(convert<Vector<std::complex<T>, 4>>(boxInfo.q2) - k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta);

    std::complex<T> tree = 0;
    for (int i = 0; i < 4; i++) {
        std::array<Vector<std::complex<T>, 4>, 3> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
        std::array<Vector<std::complex<T>, 4>, 3> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
        tree += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
            * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
    }
    counter *= tree;
    return (numerator / denominator) - counter;
}

template <typename T>
std::complex<T> gluon4vertexLoopDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta) * (square(k) + std::complex<T>(1.0i) * Delta);
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = 9.0 * square(k);
    return numerator / denominator - numerator / uvDenominator;
}

template <typename T>
std::complex<T> gluonFermionLoopMassiveDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k + q) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> denominator = (square(k) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta) * (square(k + q) - boxInfo.mass * boxInfo.mass + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = 2.0 * dot(k, k + q) - 2.0 * (dot(q, k) * dot(q, k + q) / square(q)) - 3.0 * (dot(k, k + q) - boxInfo.mass * boxInfo.mass);
    return numerator / denominator - numerator / uvDenominator;
}

template <typename T>
std::complex<T> gluonFermionLoopMasslessDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k + q) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta) * (square(k + q) + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = 2.0 * dot(k, k + q) - 2.0 * (dot(q, k) * dot(q, k + q) / square(q)) - 3.0 * (dot(k, k + q));
    return numerator / denominator - numerator / uvDenominator;
}

template <typename T>
std::complex<T> gluon3vertexLoopDiagram(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k + q) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta) * (square(k + q) + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = 3.0 * dot(q-k, k-q) + dot(q-k, -std::complex<double>(2.0)*k-q) - dot(q, q-k) * dot(q, -std::complex<double>(2.0)*k-q) / square(q) + dot(q-k, k+std::complex<double>(2.0)*q) - dot(q, q-k) * dot(q, std::complex<double>(2.0)*q+k) / square(q)
        + dot(k-q, q+std::complex<double>(2.0)*k) - dot(q, k-q) * dot(q, std::complex<double>(2.0)*k+q) / square(q) + 4.0 * dot(q + std::complex<double>(2.0)*k, -std::complex<double>(2.0)*k-q) - 4.0 / square(q) * dot(q, q+std::complex<double>(2.0)*k) * dot(q, -std::complex<double>(2.0)*k-q) + dot(q+std::complex<double>(2.0)*k, k+std::complex<double>(2.0)*q) - dot(q, q+std::complex<double>(2.0)*k) * dot(q, std::complex<double>(2.0)*q+k) / square(q)
        + dot(-k-std::complex<double>(2.0)*q, k-q) - dot(q, -std::complex<double>(2.0)*q-k) * dot(q, k-q) / square(q) + dot(-k-std::complex<double>(2.0)*q, -std::complex<double>(2.0)*k-q) - dot(q, -std::complex<double>(2.0)*q-k) * dot(q, -std::complex<double>(2.0)*k-q) / square(q) + 3.0 * dot(-k-std::complex<double>(2.0)*q, k+std::complex<double>(2.0)*q);
    return numerator / denominator - numerator / uvDenominator;
}

template <typename T>
std::complex<T> gluonFermionLoopMasslessDiagramV2(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta) * (square(k + q) + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = square(q);
    return numerator / denominator - numerator / uvDenominator;
}

template <typename T>
std::complex<T> gluonFermionLoopMassiveDiagramV2(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta) * (square(k + q) + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = square(q) + boxInfo.mass * boxInfo.mass * (2.0 + boxInfo.mass * boxInfo.mass / square(q));
    return numerator / denominator - numerator / uvDenominator;
}

template <typename T>
std::complex<T> gluonLoopCAV2(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta) * (square(k + q) + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = square(q);// * (1.0 / 8.0 - 25.0 / 64.0 + 3.0);
    return numerator / denominator - numerator / uvDenominator;
}

template <typename T>
std::complex<T> oneLoopAmplitude(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    T t1 = boxInfo.t1;
    T t2 = boxInfo.t2;
    T ca = boxInfo.ca;
    T cf = boxInfo.cf;
    T piFactor = 1.0 / ::pow(2.0 * M_PI, 4);
    std::complex<T> tree = 0;
    for (int i = 0; i < 4; i++) {
        std::array<Vector<std::complex<T>, 4>, 3> momentaMassless = { convert<Vector<std::complex<T>, 4>>(gammaUp[i]) };
        std::array<Vector<std::complex<T>, 4>, 3> momentaMassive = { convert<Vector<std::complex<T>, 4>>(gammaDown[i]) };
        tree += computeSpinorProductMassless<T>(true, boxInfo.p2, boxInfo.p1, boxInfo.spinP2, boxInfo.spinP1, momentaMassless, 1)
            * computeSpinorProductMassless<T>(false, boxInfo.q1, boxInfo.q2, boxInfo.spinQ1, boxInfo.spinQ2, momentaMassive, 1);
    }
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    T s = square(q).real();
    std::complex<T> result = -1.0 / (4.0 * s) * (ca - 2.0 * cf) * piFactor * (t1 - 1.0 / ca * t2) * masslessTriangleDiagram<T>(v, boxInfo, Delta, direction);
    result += 1.0 / (4.0 * s) * (ca - 2.0 * cf) * piFactor * (t1 - 1.0 / ca * t2) * massiveTriangleDiagram<T>(v, boxInfo, Delta, direction);
    result += 1.0 / (4.0 * s) * (ca) * piFactor * (t1 - 1.0 / ca * t2) * masslessGluonTriangleDiagram<T>(v, boxInfo, Delta, direction);
    result += 1.0 / (4.0 * s) * (ca) * piFactor * (t1 - 1.0 / ca * t2) * massiveGluonTriangleDiagram<T>(v, boxInfo, Delta, direction);
    result += (-1.0 / (4.0 * s) * (ca - 2.0 * cf) * (t1 - 1.0 / ca * t2) + 1.0 / (2.0 * s) * t1) * piFactor * boxDiagram<T>(v, boxInfo, Delta, direction);
    result += (-1.0 / (4.0 * s) * (ca - 2.0 * cf) * (t1 - 1.0 / ca * t2) - 1.0 / (4.0 * s) * ((ca - 2.0 * cf) * t1 - t2)) * piFactor * crossBoxDiagram<T>(v, boxInfo, Delta, direction);
    std::complex<T> gluonSelf = -2.0 / 3.0 / s / s * cf * (t1 - 1.0 / ca * t2) * piFactor * tree * gluonFermionLoopMasslessDiagramV2<T>(v, boxInfo, Delta, direction)
        += -2.0 / 3.0 / s / s * cf * (t1 - 1.0 / ca * t2) * piFactor * tree * gluonFermionLoopMassiveDiagramV2<T>(v, boxInfo, Delta, direction)
        += -1.0 / 12.0 * (-5.0) / s / s * ca * (t1 - 1.0 / ca * t2) * piFactor * tree * gluonLoopCAV2<T>(v, boxInfo, Delta, direction);
    result += gluonSelf;
    /*std::complex<T> gluonSelf = -1.0 / 2.0 * ca * gluon4vertexLoopDiagram<T>(v, boxInfo, Delta, direction);
    gluonSelf += -ca * gluon3vertexLoopDiagram<T>(v, boxInfo, Delta, direction);
    gluonSelf += -4 * cf * gluonFermionLoopMasslessDiagram<T>(v, boxInfo, Delta, direction);
    gluonSelf += -4 * cf * gluonFermionLoopMassiveDiagram<T>(v, boxInfo, Delta, direction);
    result += 1.0 / (6.0 * s * s) * piFactor * tree * (t1 - 1.0 / ca * t2) * gluonSelf;*/
    //std::cout << result << std::endl;
    return result * std::complex<T>(4 * 3.14 * 3.14);
}

template <typename T>
std::complex<T> cutoffTest(Vector<T, 4> v, ProcessInfo<T> boxInfo) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * std::complex<T>(1.0i), v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    //return ::exp(-v[0]*v[0] -v[1] * v[1] - v[2] * v[2] - v[3] * v[3]);
    std::complex<T> delta = 0.0i;
    std::complex<T> term = 1.0 / (square(k) + delta) / (square(q - k) + delta);
    std::complex<T> counter = 1.0 / (square(k) - boxInfo.uvScale * boxInfo.uvScale + delta) / (square(k) - boxInfo.uvScale * boxInfo.uvScale + delta);
    //if (v[0] > 0.0 && v[0] < 1.0)// && (v[0] > 0.1 || v[0] < -0.1))
        return std::complex<T>(1.0i) * (term - counter);
    //else
        //return 0.0;
}

template <std::complex<double> F(Vector<std::complex<double>, 4>, ProcessInfo<double>, double, std::complex<double>)>
int integratePolesDouble(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale = sphericalMapScale3D(x);
    Vector<double, 4> kReal = sphericalMap3D<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]), kReal[1], kReal[2], kReal[3]};
    double delta = 0.1;
    std::complex<double> direction = std::complex<double>(1.0);

    Vector<std::complex<double>, 4> q = convert<Vector<std::complex<double>, 4>>(static_cast<ProcessInfo<double>*>(userdata)->p1) + convert<Vector<std::complex<double>, 4>>(static_cast<ProcessInfo<double>*>(userdata)->p2);
    // for any diagram, the rectangle must be at least q0 wide to include all required poles
    double q0 = q[0].real();
    double rectDelta = 0.1;
    double length = 2.0 * q0 + 4.0 * rectDelta;
    double y = length * kReal[0];
    std::complex<double> dx;
    if (y < 2.0 * rectDelta) {
       k[0] = y * std::complex<double>(1.0i);
       dx = std::complex<double>(1.0i) * length;
    } else if (y < 2.0 * rectDelta + q0) {
       k[0] = std::complex<double>(2.0i) * rectDelta + (y - 2.0 * rectDelta);
       dx = length;
    } else if (y < 4.0 * rectDelta + q0) {
       k[0] = std::complex<double>(2.0i) * rectDelta + q0 - std::complex<double>(1.0i) * (y - 2.0 * rectDelta - q0);
       dx = -std::complex<double>(1.0i) * length;
    } else {
       k[0] = q0 - (y - 4.0 * rectDelta - q0);
       dx = -1.0 * length;
    }

    // minus since contour is wrong way around
    std::complex<double> result = -dx * F(k, *static_cast<ProcessInfo<double>*>(userdata), delta, direction) * scale;

    if (y < 2.0 * rectDelta) {
       k[0] = -std::complex<double>(2.0i) * rectDelta - q0 + y * std::complex<double>(1.0i);
       dx = std::complex<double>(1.0i) * length;
    } else if (y < 2.0 * rectDelta + q0) {
       k[0] = -std::complex<double>(2.0i) * rectDelta - q0 + std::complex<double>(2.0i) * rectDelta + (y - 2.0 * rectDelta);
       dx = length;
    } else if (y < 4.0 * rectDelta + q0) {
       k[0] = -std::complex<double>(2.0i) * rectDelta - q0 + std::complex<double>(2.0i) * rectDelta + q0 - std::complex<double>(1.0i) * (y - 2.0 * rectDelta - q0);
       dx = -std::complex<double>(1.0i) * length;
    } else {
       k[0] = -std::complex<double>(2.0i) * rectDelta - q0 + q0 - (y - 4.0 * rectDelta - q0);
       dx = -1.0 * length;
    }
    result += dx * F(k, *static_cast<ProcessInfo<double>*>(userdata), delta, direction) * scale;

    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

template <std::complex<double> F(Vector<double, 4>, ProcessInfo<double>)>
int integrateDouble3D(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale = sphericalMapScale3D(x);
    Vector<double, 4> k = sphericalMap3D<double>(const_cast<double*>(x));
    std::complex<double> result = F(k, *static_cast<ProcessInfo<double>*>(userdata)) * scale;
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

template <typename T>
std::complex<T> residue(Vector<T, 4> v, ProcessInfo<T> boxInfo) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]), v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    T sRoot = ::sqrt(square(q).real());
    T delta = 0.1;
    T length = 2.0 * sRoot + 4.0 * delta;
    T x = length * v[0];
    std::complex<T> dx;
    if (x < 2.0 * delta) {
       k[0] = x * std::complex<T>(1.0i);
       dx = std::complex<T>(1.0i) * length;
    } else if (x < 2.0 * delta + sRoot) {
       k[0] = std::complex<T>(2.0i) * delta + (x - 2.0 * delta);
       dx = length;
    } else if (x < 4.0 * delta + sRoot) {
       k[0] = std::complex<T>(2.0i) * delta + sRoot - std::complex<T>(1.0i) * (x - 2.0 * delta - sRoot);
       dx = -std::complex<T>(1.0i) * length;
    } else {
       k[0] = sRoot - (x - 4.0 * delta - sRoot);
       dx = -1.0 * length;
    }
    //k[0] = k[0] + std::complex<T>(0.00001i);
    //k[0] = k[0] - std::complex<T>(1.0i) * delta;

    //std::complex<T> result = ::exp(-v[1] * v[1] - v[2] * v[2] - v[3] * v[3]) * k[0] * dx;
    //return result;

    std::complex<T> term = 1.0 / (square(k) + std::complex<T>(0.0001i)/*+ std::complex<T>(1.0i) * delta*/) / (square(q - k) + std::complex<T>(0.0001i) /*+ std::complex<T>(1.0i) * delta*/);
    std::complex<T> counter = 1.0 / (square(k) - boxInfo.uvScale * boxInfo.uvScale /*+ std::complex<T>(1.0i) * delta*/) / (square(k) - boxInfo.uvScale * boxInfo.uvScale /*+ std::complex<T>(1.0i) * delta*/);
    //if (v[0] > 0.0 && v[0] < 1.0)// && (v[0] > 0.1 || v[0] < -0.1))
    // minus since contour is wrong way around
        return -(term - counter) * dx;
    //else
        //return 0.0;
}

template <typename T>
std::complex<T> residue3D(Vector<T, 4> v, ProcessInfo<T> boxInfo) {
    std::complex<T> k2 = v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    T s = square(q).real();
    std::complex<T> delta = 0.0001i;
    if (s <= k2.real()) {
       return 0;
    } else {
        auto constexpr pi = M_PI;
        return std::complex<T>(2.0i) * pi / (s - 2.0 * ::sqrt(s*k2.real()) + delta + delta * 2.0 * (::sqrt(s) - ::sqrt(k2.real()))) / (-std::complex<T>(2.0) * ::sqrt(k2.real()) + delta * 2.0);
    }
}

template <typename T>
std::complex<T> test(Vector<std::complex<T>, 4> v, ProcessInfo<T> boxInfo, double Delta, std::complex<double> direction) {
    Vector<std::complex<T>, 4> k = {std::complex<T>(v[0]) * direction, v[1], v[2], v[3]};
    Vector<std::complex<T>, 4> q = convert<Vector<std::complex<T>, 4>>(boxInfo.p1) + convert<Vector<std::complex<T>, 4>>(boxInfo.p2);
    //std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    //std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta) * (square(k + q) + std::complex<T>(1.0i) * Delta);
    std::complex<T> uvDenominator = (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta) * (square(k) - boxInfo.uvScale * boxInfo.uvScale + std::complex<T>(1.0i) * Delta);
    std::complex<T> denominator = (square(k) + std::complex<T>(1.0i) * Delta) * (square(k + q) + std::complex<T>(1.0i) * Delta);
    std::complex<T> numerator = 2.0 * dot(k, q);
    //std::complex<T> numerator = 2.0 * dot(k, k + q) - 2.0 * (dot(q, k) * dot(q, k + q) / square(q)) - 3.0 * (dot(k, k + q));
    return numerator / denominator - numerator / uvDenominator;
}

int doubleLoop(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale(x);
    double scale2 = sphericalMapScale(&x[4]);
    Vector<double, 4> kReal = sphericalMap<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<double, 4> lReal = sphericalMap<double>(const_cast<double*>(&x[4]));
    Vector<std::complex<double>, 4> l = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> p = {std::complex<double>(10.0), 0.0, 0.0, 0.0};
    double scale = scale1 * scale2;
    std::complex<double> delta = 0.0000001i;
    std::complex<double> integrand = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k - l) + delta);
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = -scale * integrand;
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

int doubleLoopResidue1(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale3D(x);
    double scale2 = sphericalMapScale(&x[4]);
    Vector<double, 4> kReal = sphericalMap3D<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<double, 4> lReal = sphericalMap<double>(const_cast<double*>(&x[4]));
    Vector<std::complex<double>, 4> l = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    //Vector<std::complex<double>, 4> l1 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> p = {std::complex<double>(10.0), 0.0, 0.0, 0.0};

    double rectDelta = 0.1;
    double squareLength = rectDelta + std::sqrt(square(p).real());
    double length = 8.0 * squareLength;
    double y = length * kReal[0];
    std::complex<double> dx;
    if (y < squareLength) {
       k[0] = y * std::complex<double>(1.0);
       dx = std::complex<double>(1.0) * length;
    } else if (y < 2.0 * squareLength) {
       k[0] = std::complex<double>(1.0) * squareLength + std::complex<double>(1.0i) * (y - 1.0 * squareLength);
       dx = std::complex<double>(1.0i) * length;
    } else if (y < 3.0 * squareLength) {
       k[0] = std::complex<double>(1.0 + 1.0i) * squareLength + std::complex<double>(-1.0) * (y - 2.0 * squareLength);
       dx = std::complex<double>(-1.0) * length;
    } else if (y < 5.0 * squareLength) {
       k[0] = std::complex<double>(1.0i) * squareLength + std::complex<double>(-1.0i) * (y - 3.0 * squareLength);
       dx = std::complex<double>(-1.0i) * length;
    } else if (y < 6.0 * squareLength) {
       k[0] = std::complex<double>(-1.0i) * squareLength + std::complex<double>(-1.0) * (y - 5.0 * squareLength);
       dx = std::complex<double>(-1.0) * length;
    } else if (y < 7.0 * squareLength) {
       k[0] = std::complex<double>(-1.0 - 1.0i) * squareLength + std::complex<double>(1.0i) * (y - 6.0 * squareLength);
       dx = std::complex<double>(1.0i) * length;
    } else  {
       k[0] = std::complex<double>(-1.0) * squareLength + std::complex<double>(1.0) * (y - 7.0 * squareLength);
       dx = std::complex<double>(1.0) * length;
    }
    //std::cout << "x " << squareLength << std::endl;
    /*
    std::cout << "k " << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << std::endl;
    std::cout << "l " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
    */
    double scale = scale1 * scale2;
    std::complex<double> delta = 0.001i;
    std::complex<double> integrand = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k - l) + delta);
    //std::complex<double> integrand1 = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l1) + delta) / (square(p - l1) + delta) / (square(k - l1) + delta);
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = dx * scale * (integrand * std::complex<double>(1.0i));
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

int doubleLoopResidueAnalytic(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale3D(x);
    double scale2 = sphericalMapScale(&x[4]);
    Vector<double, 4> kReal = sphericalMap3D<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<double, 4> lReal = sphericalMap<double>(const_cast<double*>(&x[4]));
    Vector<std::complex<double>, 4> l = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    //Vector<std::complex<double>, 4> l1 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> p = {std::complex<double>(10.0), 0.0, 0.0, 0.0};

    double scale = scale1 * scale2;
    std::complex<double> delta = 0.1i;
    k[0] = p[0] - std::sqrt(kReal[1] * kReal[1] * kReal[2] * kReal[2] * kReal[3] * kReal[3]) + delta;
    std::complex<double> integrand = std::complex<double>(2.0 * 3.14i) / (square(k) + delta) / (-2.0 * std::sqrt(kReal[1] * kReal[1] * kReal[2] * kReal[2] * kReal[3] * kReal[3]) + 2.0 * delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k - l) + delta);
    //std::complex<double> integrand1 = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l1) + delta) / (square(p - l1) + delta) / (square(k - l1) + delta);
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = scale * (integrand * std::complex<double>(1.0i));
    if (p[0].real() <= std::sqrt(kReal[1] * kReal[1] * kReal[2] * kReal[2] * kReal[3] * kReal[3])) {
        result = 0;
    }
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

int doubleLoopResidueAnalytic2(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale3D(x);
    double scale2 = sphericalMapScale3D(&x[4]);
    Vector<double, 4> kReal = sphericalMap3D<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<double, 4> lReal = sphericalMap3D<double>(const_cast<double*>(&x[4]));
    Vector<std::complex<double>, 4> l = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    //Vector<std::complex<double>, 4> l1 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> p = {std::complex<double>(10.0), 0.0, 0.0, 0.0};

    double scale = scale1 * scale2;
    std::complex<double> delta = 0.001i;
    k[0] = p[0] - std::sqrt(kReal[1] * kReal[1] * kReal[2] * kReal[2] * kReal[3] * kReal[3]) + delta;
    l[0] = p[0] - std::sqrt(lReal[1] * lReal[1] * lReal[2] * lReal[2] * lReal[3] * lReal[3]) + delta;
    std::complex<double> integrand = std::complex<double>(2.0 * 3.14i) * std::complex<double>(2.0 * 3.14i) / (square(k) + delta) / (-2.0 * std::sqrt(kReal[1] * kReal[1] * kReal[2] * kReal[2] * kReal[3] * kReal[3]) + 2.0 * delta) / (square(l) + delta) / (-2.0 * std::sqrt(lReal[1] * lReal[1] * lReal[2] * lReal[2] * lReal[3] * lReal[3]) + 2.0 * delta) / (square(k - l) + delta);
    //std::complex<double> integrand1 = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l1) + delta) / (square(p - l1) + delta) / (square(k - l1) + delta);
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = scale * (integrand);
    if (p[0].real() <= std::sqrt(kReal[1] * kReal[1] * kReal[2] * kReal[2] * kReal[3] * kReal[3])) {
        result = 0;
    }
    if (p[0].real() <= std::sqrt(lReal[1] * lReal[1] * lReal[2] * lReal[2] * lReal[3] * lReal[3])) {
        result = 0;
    }
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

int doubleLoopResidue1v2(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale3D(x);
    double scale2 = sphericalMapScale(&x[4]);
    Vector<double, 4> kReal = sphericalMap3D<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k1 = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> k2 = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> k3 = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> k4 = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<double, 4> lReal = sphericalMap<double>(const_cast<double*>(&x[4]));
    Vector<std::complex<double>, 4> l = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    //Vector<std::complex<double>, 4> l1 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> p = {std::complex<double>(10.0), 0.0, 0.0, 0.0};

    double rectDelta = 0.1;
    double squareLength = rectDelta + std::sqrt(square(p).real());
    double length = squareLength;
    double y = length * kReal[0];
    std::complex<double> dx1;
    std::complex<double> dx2;
    std::complex<double> dx3;
    std::complex<double> dx4;
    std::complex<double> delta = 0.0000001i;
    k1[0] = y;
    k2[0] = delta * 2.0 + y;
    dx1 = 1.0;
    dx2 = -1.0;

    k3[0] = -y;
    k4[0] = -delta * 2.0 - y;
    dx3 = -1.0;
    dx4 = 1.0;
    //std::cout << "x " << squareLength << std::endl;
    /*
    std::cout << "k " << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << std::endl;
    std::cout << "l " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
    */
    double scale = scale1 * scale2;
    std::complex<double> integrand1 = dx1 * 1.0 / (square(k1) + delta) / (square(p - k1) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k1 - l) + delta);
    std::complex<double> integrand2 = dx2 * 1.0 / (square(k2) + delta) / (square(p - k2) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k2 - l) + delta);
    std::complex<double> integrand3 = dx3 * 1.0 / (square(k3) + delta) / (square(p - k3) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k3 - l) + delta);
    std::complex<double> integrand4 = dx3 * 1.0 / (square(k4) + delta) / (square(p - k4) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k4 - l) + delta);
    //std::complex<double> integrand1 = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l1) + delta) / (square(p - l1) + delta) / (square(k - l1) + delta);
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = scale * (integrand1+integrand2+integrand3+integrand4) * std::complex<double>(1.0i);
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

std::complex<double> doubleLoopIntegrandv2(Vector<std::complex<double>, 4> k, Vector<std::complex<double>, 4> l) {
    Vector<std::complex<double>, 4> p = {std::complex<double>(10.0), 0.0, 0.0, 0.0};
    std::complex<double> delta = 0.0000001i;
    std::complex<double> integrand = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k - l) + delta);
    return integrand;
}

int doubleLoopResidue2v2(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale3D(x);
    double scale2 = sphericalMapScale3D(&x[4]);
    Vector<double, 4> kReal = sphericalMap3D<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k1 = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> k2 = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> k3 = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> k4 = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<double, 4> lReal = sphericalMap3D<double>(const_cast<double*>(&x[4]));
    Vector<std::complex<double>, 4> l1 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> l2 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> l3 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> l4 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    //Vector<std::complex<double>, 4> l1 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> p = {std::complex<double>(10.0), 0.0, 0.0, 0.0};

    double rectDelta = 0.001;
    double squareLength = rectDelta + std::sqrt(square(p).real());
    double length = squareLength;

    //std::cout << "x " << squareLength << std::endl;
    /*
    std::cout << "k " << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << std::endl;
    std::cout << "l " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
    */
    double scale = scale1 * scale2;
    std::complex<double> delta = 0.0000001i;


    double y = length * kReal[0];
    std::complex<double> dk1;
    std::complex<double> dk2;
    std::complex<double> dk3;
    std::complex<double> dk4;
    k1[0] = y;
    k2[0] = delta * 2.0 + y;
    dk1 = 1.0;
    dk2 = -1.0;
    k3[0] = -y;
    k4[0] = -delta * 2.0 - y;
    dk3 = -1.0;
    dk4 = 1.0;


    y = length * lReal[0];
    std::complex<double> dl1;
    std::complex<double> dl2;
    std::complex<double> dl3;
    std::complex<double> dl4;
    l1[0] = y;
    l2[0] = delta * 2.0 + y;
    dl1 = 1.0;
    dl2 = -1.0;
    l3[0] = -y;
    l4[0] = -delta * 2.0 - y;
    dl3 = -1.0;
    dl4 = 1.0;



    std::complex<double> integrand1 = dk1 * dl1 * doubleLoopIntegrandv2(k1, l1)
        + dk1 * dl2 * doubleLoopIntegrandv2(k1, l2);
        + dk1 * dl3 * doubleLoopIntegrandv2(k1, l3);
        + dk1 * dl4 * doubleLoopIntegrandv2(k1, l4);
    std::complex<double> integrand2 = dk2 * dl1 * doubleLoopIntegrandv2(k2, l1)
        + dk2 * dl2 * doubleLoopIntegrandv2(k2, l2);
        + dk2 * dl3 * doubleLoopIntegrandv2(k2, l3);
        + dk2 * dl4 * doubleLoopIntegrandv2(k2, l4);
    std::complex<double> integrand3 = dk3 * dl1 * doubleLoopIntegrandv2(k3, l1)
        + dk3 * dl2 * doubleLoopIntegrandv2(k3, l2);
        + dk3 * dl3 * doubleLoopIntegrandv2(k3, l3);
        + dk3 * dl4 * doubleLoopIntegrandv2(k3, l4);
    std::complex<double> integrand4 = dk4 * dl1 * doubleLoopIntegrandv2(k4, l1)
        + dk4 * dl2 * doubleLoopIntegrandv2(k4, l2);
        + dk4 * dl3 * doubleLoopIntegrandv2(k4, l3);
        + dk4 * dl4 * doubleLoopIntegrandv2(k4, l4);

    //std::complex<double> integrand1 = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l1) + delta) / (square(p - l1) + delta) / (square(k - l1) + delta);
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = scale * (integrand1 + integrand2 + integrand3 + integrand4);
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

int doubleLoopResidue2(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale3D(x);
    double scale2 = sphericalMapScale3D(&x[4]);
    Vector<double, 4> kReal = sphericalMap3D<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> k2 = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<double, 4> lReal = sphericalMap3D<double>(const_cast<double*>(&x[4]));
    Vector<std::complex<double>, 4> l = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> l2 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    //Vector<std::complex<double>, 4> l1 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> p = {std::complex<double>(1.0), 0.0, 0.0, 0.0};

    double rectDelta = 0.1;
    double squareLength = rectDelta + std::sqrt(square(p).real());
    double length = 4.0 * squareLength;


    double a = length * kReal[0];
    std::complex<double> da1;
    std::complex<double> da2;
    if (a < squareLength) {
       k[0] = a * std::complex<double>(1.0);
       da1 = std::complex<double>(1.0) * length;
       k2[0] = a * std::complex<double>(-1.0);
       da2 = std::complex<double>(1.0) * length;
    } else if (a < 2.0 * squareLength) {
       k[0] = std::complex<double>(1.0) * squareLength + std::complex<double>(1.0i) * (a - 1.0 * squareLength);
       da1 = std::complex<double>(1.0i) * length;
       k2[0] = std::complex<double>(-1.0) * squareLength + std::complex<double>(-1.0i) * (a - 1.0 * squareLength);
       da2 = std::complex<double>(1.0i) * length;
    } else if (a < 3.0 * squareLength) {
       k[0] = std::complex<double>(1.0 + 1.0i) * squareLength + std::complex<double>(-1.0) * (a - 2.0 * squareLength);
       da1 = std::complex<double>(-1.0) * length;
       k2[0] = std::complex<double>(-1.0 - 1.0i) * squareLength + std::complex<double>(1.0) * (a - 2.0 * squareLength);
       da2 = std::complex<double>(-1.0) * length;
    } else if (a < 4.0 * squareLength) {
       k[0] = std::complex<double>(1.0i) * squareLength + std::complex<double>(-1.0i) * (a - 3.0 * squareLength);
       da1 = std::complex<double>(-1.0i) * length;
       k2[0] = std::complex<double>(-1.0i) * squareLength + std::complex<double>(1.0i) * (a - 3.0 * squareLength);
       da2 = std::complex<double>(-1.0i) * length;
    }

    double b = length * kReal[0];
    std::complex<double> db1;
    std::complex<double> db2;
    if (b < squareLength) {
       l[0] = b * std::complex<double>(1.0);
       db1 = std::complex<double>(1.0) * length;
       l2[0] = b * std::complex<double>(-1.0);
       db2 = std::complex<double>(1.0) * length;
    } else if (b < 2.0 * squareLength) {
       l[0] = std::complex<double>(1.0) * squareLength + std::complex<double>(1.0i) * (b - 1.0 * squareLength);
       db1 = std::complex<double>(1.0i) * length;
       l2[0] = std::complex<double>(-1.0) * squareLength + std::complex<double>(-1.0i) * (b - 1.0 * squareLength);
       db2 = std::complex<double>(1.0i) * length;
    } else if (b < 3.0 * squareLength) {
       l[0] = std::complex<double>(1.0 + 1.0i) * squareLength + std::complex<double>(-1.0) * (b - 2.0 * squareLength);
       db1 = std::complex<double>(-1.0) * length;
       l2[0] = std::complex<double>(-1.0 - 1.0i) * squareLength + std::complex<double>(1.0) * (b - 2.0 * squareLength);
       db2 = std::complex<double>(-1.0) * length;
    } else if (b < 4.0 * squareLength) {
       l[0] = std::complex<double>(1.0i) * squareLength + std::complex<double>(-1.0i) * (b - 3.0 * squareLength);
       db1 = std::complex<double>(-1.0i) * length;
       l2[0] = std::complex<double>(-1.0i) * squareLength + std::complex<double>(1.0i) * (b - 3.0 * squareLength);
       db2 = std::complex<double>(-1.0i) * length;
    }
    //std::cout << "x " << squareLength << std::endl;
    /*
    std::cout << "k " << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << std::endl;
    std::cout << "l " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
    */
    double scale = scale1 * scale2;
    std::complex<double> delta = 0.001i;
    std::complex<double> integrand1 = da1 * db1 * 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k - l) + delta);
    std::complex<double> integrand2 = da2 * db2 * 1.0 / (square(k2) + delta) / (square(p - k2) + delta) / (square(l2) + delta) / (square(p - l2) + delta) / (square(k2 - l2) + delta);
    std::complex<double> integrand3 = da1 * db2 * 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l2) + delta) / (square(p - l2) + delta) / (square(k - l2) + delta);
    std::complex<double> integrand4 = da2 * db1 * 1.0 / (square(k2) + delta) / (square(p - k2) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k2 - l) + delta);
    //std::complex<double> integrand1 = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l1) + delta) / (square(p - l1) + delta) / (square(k - l1) + delta);
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = scale * (integrand1 + integrand2 + integrand3 + integrand4);
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}










std::complex<double> doubleLoopIntegrand(double scale, Vector<std::complex<double>, 4> k, Vector<std::complex<double>, 4> l, Vector<std::complex<double>, 4> l1) {
    Vector<std::complex<double>, 4> p = {std::complex<double>(1.0), 0.0, 0.0, 0.0};
    std::complex<double> delta = 0.01i * scale;
    std::complex<double> integrand = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l) + delta) / (square(p - l) + delta) / (square(k - l) + delta);
    std::complex<double> integrand1 = 1.0 / (square(k) + delta) / (square(p - k) + delta) / (square(l1) + delta) / (square(p - l1) + delta) / (square(k - l1) + delta);
    return integrand + std::complex<double>(1.0i) * integrand1;
}

int doubleLoopResidue2bla(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale3D(x);
    double scale2 = sphericalMapScale(&x[4]);
    Vector<double, 4> kReal = sphericalMap3D<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]), kReal[1], kReal[2], kReal[3]};
    Vector<double, 4> lReal = sphericalMap<double>(const_cast<double*>(&x[4]));
    Vector<std::complex<double>, 4> l = {std::complex<double>(lReal[0]), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> l1 = {std::complex<double>(lReal[0]) * std::complex<double>(1.0i), lReal[1], lReal[2], lReal[3]};

    Vector<std::complex<double>, 4> p = {std::complex<double>(1.0), 0.0, 0.0, 0.0};
    double rectDelta = 0.1;
    double squareLength = rectDelta + std::sqrt(square(p).real()) + std::abs(l[0]) + std::sqrt(std::abs(l[0] - square(l) + k[1] * k[1] + k[2] * k[2] + k[3] * k[3] - 2.0 * k[1] * l[1] + k[2] * l[2] + k[3] * l[3]));
    //double squareLength = rectDelta + 2.0 * std::abs(l[0]);
    //squareLength = 1.0;
    //std::cout << l[0] << std::endl;
    double length = squareLength;
    double y = squareLength * kReal[0];
    double z = kReal[0];
    std::complex<double> sum;

    double scale = scale1 * scale2;
    std::complex<double> dx;
       k[0] = y * std::complex<double>(1.0);
       dx = std::complex<double>(1.0) * length;
       sum += dx * doubleLoopIntegrand(scale, k, l, l1);

       /*k[0] = std::complex<double>(1.0) * squareLength + std::complex<double>(1.0i) * (y);
       dx = std::complex<double>(1.0i) * length;
       sum += dx * doubleLoopIntegrand(k, l, l1);
       */

       k[0] = std::complex<double>(1.0 + 1.0i) * squareLength + std::complex<double>(-1.0) * (squareLength - y);
       dx = std::complex<double>(-1.0) * length;
       sum += dx * doubleLoopIntegrand(scale, k, l, l1);

       /*k[0] = std::complex<double>(1.0i) * squareLength + std::complex<double>(-1.0i) * (squareLength - y);
       dx = std::complex<double>(-1.0i) * length;
       sum += dx * doubleLoopIntegrand(k, l, l1);
       */
//  ------------------------------------------------------------------------------------------------------------
       /*
       k[0] = std::complex<double>(1.0i) * squareLength + std::complex<double>(-1.0i) * (y);
       dx = std::complex<double>(-1.0i) * length;
       sum += dx * doubleLoopIntegrand(k, l, l1);

       k[0] = std::complex<double>(-1.0i) * squareLength + std::complex<double>(-1.0) * (y);
       dx = std::complex<double>(-1.0) * length;
       sum += dx * doubleLoopIntegrand(k, l, l1);

       k[0] = std::complex<double>(-1.0 - 1.0i) * squareLength + std::complex<double>(1.0i) * (squareLength - y);
       dx = std::complex<double>(1.0i) * length;
       sum += dx * doubleLoopIntegrand(k, l, l1);

       k[0] = std::complex<double>(-1.0) * squareLength + std::complex<double>(1.0) * (squareLength - y);
       dx = std::complex<double>(1.0) * length;
       sum += dx * doubleLoopIntegrand(k, l, l1);
       */
    //std::cout << "x " << squareLength << std::endl;
    /*
    std::cout << "k " << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << std::endl;
    std::cout << "l " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
    */
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = dx * scale * (sum);
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

int doubleLoop2(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale(x);
    double scale2 = sphericalMapScale(&x[4]);
    double scale = scale1;// * scale2;
    Vector<double, 4> kReal = sphericalMap<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> k2 = {std::complex<double>(kReal[0]), kReal[1], kReal[2], kReal[3]};
    bool fail = false;
    if (std::abs(square(k) - 1.0) < 0.1 * scale) {
        fail = true;
        std::complex<double> length = k[1] * k[1] * k[2] * k[2] * k[3] * k[3];
        k2[0] = k[0] / std::abs(k[0]) * length;
        k2[1] = k[1] / length * std::abs(k[0]);
        k2[2] = k[2] / length * std::abs(k[0]);
        k2[3] = k[3] / length * std::abs(k[0]);
    }
    Vector<double, 4> lReal = sphericalMap<double>(const_cast<double*>(&x[4]));
    Vector<std::complex<double>, 4> l = {std::complex<double>(lReal[0]), lReal[1], lReal[2], lReal[3]};
    Vector<std::complex<double>, 4> l2 = {std::complex<double>(lReal[0]), lReal[1], lReal[2], lReal[3]};
    if (std::abs(square(l) - 1.0) < 0.1 * scale) {
        //fail = true;
        std::complex<double> length = l[1] * l[1] * l[2] * l[2] * l[3] * l[3];
        l2[0] = l[0] / std::abs(l[0]) * length;
        l2[1] = l[1] / length * std::abs(l[0]);
        l2[2] = l[2] / length * std::abs(l[0]);
        l2[3] = l[3] / length * std::abs(l[0]);
    }
    Vector<std::complex<double>, 4> p = {std::complex<double>(1.0), 0.0, 0.0, 0.0};
    std::complex<double> delta = 0.00i;
    std::complex<double> integrand = 1.0 / (square(k) - 1.4 + delta) / (square(k) - 1.3 + delta) / (square(k) - 1.1 + delta) / (square(k) - 1.2 + delta);// / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta);
    //integrand += 1.0 / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta);
    //integrand += 1.0 / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta) / (square(l2) - 1.0 + delta) / (square(l2) - 1.0 + delta) / (square(l2) - 1.0 + delta);
    //integrand += 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l2) - 1.0 + delta) / (square(l2) - 1.0 + delta) / (square(l2) - 1.0 + delta);
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = scale * integrand;
    if (fail) {
        //std::cout << scale << std::endl;
        //result = 0.0;
    } else {
        //std::cout << scale << std::endl;
    }
        //std::cout << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << " " << std::endl;
    //std::cout << integrand << std::endl;
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

int doubleLoop24d(const int* ndim, const double x[], const int* ncomp, double f[], void* userdata) {
    double scale1 = sphericalMapScale(x);
    double scale = scale1;// * scale2;
    Vector<double, 4> kReal = sphericalMap<double>(const_cast<double*>(x));
    Vector<std::complex<double>, 4> k = {std::complex<double>(kReal[0]) * std::complex<double>(1.0i), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> k2 = {std::complex<double>(kReal[0]), kReal[1], kReal[2], kReal[3]};
    Vector<std::complex<double>, 4> p = {std::complex<double>(1.0), 0.0, 0.0, 0.0};
        //std::cout << k[0] << " " << k[1] << " " << k[2] << " " << k[3] << " " << std::endl;
        //std::cout << scale << std::endl;
    std::complex<double> delta = 0.1i;
    /*if (std::abs(square(k) - 1.0) < 0.01) {
        return 0;
        std::complex<double> length = k[1] * k[1] * k[2] * k[2] * k[3] * k[3];
        k2[0] = k[0] / std::abs(k[0]) * length;
        k2[1] = k[1] / length * std::abs(k[0]);
        k2[2] = k[2] / length * std::abs(k[0]);
        k2[3] = k[3] / length * std::abs(k[0]);
    }*/
    std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta);// / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta);
    //integrand += 1.0 / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta);// / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta);
    //std::cout << square(k) - 1.0 << std::endl;
    //std::complex<double> integrand = 1.0 / std::exp(kReal[0] * kReal[0] + kReal[1] * kReal[1] + kReal[2] * kReal[2] + kReal[3] * kReal[3]);
    //integrand += 1.0 / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta);
    //integrand += 1.0 / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta) / (square(k2) - 1.0 + delta) / (square(l2) - 1.0 + delta) / (square(l2) - 1.0 + delta) / (square(l2) - 1.0 + delta);
    //integrand += 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l2) - 1.0 + delta) / (square(l2) - 1.0 + delta) / (square(l2) - 1.0 + delta);
    //std::complex<double> integrand = 1.0 / (square(k) - 1.0 + delta) / (square(k) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(l) - 1.0 + delta) / (square(k - l) + delta);
    std::complex<double> result = scale * integrand;
    f[0] = result.real();
    f[1] = result.imag();
    return 0;
}

int main() {
    ProcessInfo<double> processInfo;
    processInfo.p1 = Vector<double, 4>(std::sqrt(5)/2.0, 0, 0, std::sqrt(5)/2.0);
    processInfo.p2 = Vector<double, 4>(std::sqrt(5)/2.0, 0, 0, -std::sqrt(5)/2.0);
    processInfo.q1 = Vector<double, 4>(std::sqrt(5)/2.0, 0, 1.0 / 2.0 * std::sqrt(19.0 / 20.0), 1.0 / 4.0 / std::sqrt(5));
    processInfo.q2 = Vector<double, 4>(std::sqrt(5)/2.0, 0, -1.0 / 2.0 * std::sqrt(19.0 / 20.0), -1.0 / 4.0 / std::sqrt(5));
    processInfo.mass = 1;
    processInfo.spinP1 = Vector<std::complex<double>, 2>(1, 0);
    processInfo.spinP2 = Vector<std::complex<double>, 2>(0, 1);
    processInfo.spinQ1 = Vector<std::complex<double>, 2>(1, 0);
    processInfo.spinQ2 = Vector<std::complex<double>, 2>(0, 1);
    processInfo.uvScale = 1;
    processInfo.t1 = 2;
    processInfo.t2 = 0;
    processInfo.ca = 3;
    processInfo.cf = 4.0/3.0;

    double s = 5.0;
    double res = 0;
    Vector<std::complex<double>, 2> spins[2] = {Vector<std::complex<double>, 2>(1, 0),Vector<std::complex<double>, 2>(0, 1)};

    for (int a = 0; a < 2; a++)
    for (int b = 0; b < 2; b++)
    for (int c = 0; c < 2; c++)
    for (int d = 0; d < 2; d++) {
        std::complex<double> tree = 0;
        for (int i = 0; i < 4; i++) {
            std::array<Vector<std::complex<double>, 4>, 3> momentaMassless = { convert<Vector<std::complex<double>, 4>>(gammaUp[i]) };
            std::array<Vector<std::complex<double>, 4>, 3> momentaMassive = { convert<Vector<std::complex<double>, 4>>(gammaDown[i]) };
            tree += computeSpinorProductMassless<double>(true, processInfo.p2, processInfo.p1, spins[a], spins[b], momentaMassless, 1)
                * computeSpinorProductMassless<double>(false, processInfo.q1, processInfo.q2, spins[c], spins[d], momentaMassive, 1);
        }
        tree /= s;
        res += std::abs(tree) * std::abs(tree);
    }
    res /= 4.0;
    std::cout << "tree: " << res << std::endl;
    //return 0;


    int neval;
    int fail;
    double integral[2];
    double error[2];
    double prob[2];

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0,100);
    /*llVegas(8, 2,
        //&integratePolesDouble<oneLoopAmplitude<double>>, static_cast<void*>(&processInfo), 1,
        &doubleLoopResidueAnalytic2, static_cast<void*>(&processInfo), 1,
        0.1, 0.1,
        3, dist(rng),
        100000, 100000000000,
        100000, 100000, 800000,
        1, "", nullptr,
        &neval, &fail,
        &integral[0], &error[0], &prob[0]);*/
    /*Vegas(8, 2,
        &doubleLoop2, static_cast<void*>(&processInfo), 1,
        0.01, 0.01,
        3, dist(rng),
        100000, 10000000000,
        1000000, 100000, 800000,
        1, "", nullptr,
        &neval, &fail,
        &integral[0], &error[0], &prob[0]);*/
    res = 0;
    for (int a = 0; a < 2; a++)
    for (int b = 0; b < 2; b++)
    for (int c = 0; c < 2; c++)
    for (int d = 0; d < 2; d++) {
        std::cout << "spin iteration " << d + 2*c + 4*b + 8*a << std::endl;
        std::complex<double> tree = 0;
        for (int i = 0; i < 4; i++) {
            std::array<Vector<std::complex<double>, 4>, 3> momentaMassless = { convert<Vector<std::complex<double>, 4>>(gammaUp[i]) };
            std::array<Vector<std::complex<double>, 4>, 3> momentaMassive = { convert<Vector<std::complex<double>, 4>>(gammaDown[i]) };
            tree += computeSpinorProductMassless<double>(true, processInfo.p2, processInfo.p1, spins[a], spins[b], momentaMassless, 1)
                * computeSpinorProductMassless<double>(false, processInfo.q1, processInfo.q2, spins[c], spins[d], momentaMassive, 1);
        }
        tree /= s;
        processInfo.spinP1 = spins[a];
        processInfo.spinP2 = spins[b];
        processInfo.spinQ1 = spins[c];
        processInfo.spinQ2 = spins[d];
        Vegas(4, 2,
            //&integratePolesDouble<oneLoopAmplitude<double>>, static_cast<void*>(&processInfo), 1,
            &integrateDiagramDouble<oneLoopAmplitude<double>>, static_cast<void*>(&processInfo), 1,
            0.1, 0.1,
            3, dist(rng),
            100000, 100000000,
            10000, 1000, 8000,
            1, "", nullptr,
            &neval, &fail,
            &integral[0], &error[0], &prob[0]);
        std::complex<double> imagAxis = integral[0] + 1.0i * integral[1];
        Vegas(4, 2,
            &integratePolesDouble<oneLoopAmplitude<double>>, static_cast<void*>(&processInfo), 1,
            //&integrateDiagramDouble<oneLoopAmplitude<double>>, static_cast<void*>(&processInfo), 1,
            0.2, 0.5,
            3, dist(rng),
            10000, 100000000,
            1000, 1000, 8000,
            1, "", nullptr,
            &neval, &fail,
            &integral[0], &error[0], &prob[0]);
        std::complex<double> residue = integral[0] + 1.0i * integral[1];
        //std::complex<double> residue = 0;
        std::cout << "result = " << imagAxis + residue;
        std::complex<double> iFactor = -1.0i;

        res += (std::conj(tree) * (iFactor * (imagAxis + residue))).real() / 4.0;
    }
    std::cout << "final result summed up: " << res << std::endl;




    /*
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0,100);
    Vegas(4, 2,
        &integrateDouble3D<residue<double>>, static_cast<void*>(&processInfo), 1,
        //&integrateDiagramDouble<cutoffTest<double>>, static_cast<void*>(&processInfo), 1,
        0.01, 0.01,
        3, dist(rng),
        10000, 10000000,
        1000, 1000, 8000,
        1, "", nullptr,
        &neval, &fail,
        &integral[0], &error[0], &prob[0]);
    std::complex<double> imagAxis = integral[0] + 1.0i * integral[1];
    Vegas(4, 2,
        //&integrateDouble3D<residue<double>>, static_cast<void*>(&processInfo), 1,
        &integrateDiagramDouble<cutoffTest<double>>, static_cast<void*>(&processInfo), 1,
        0.01, 0.01,
        3, dist(rng),
        10000, 10000000,
        1000, 1000, 8000,
        1, "", nullptr,
        &neval, &fail,
        &integral[0], &error[0], &prob[0]);
    std::complex<double> residue = integral[0] + 1.0i * integral[1];
    std::cout << "result = " << imagAxis + residue;
    */
}
