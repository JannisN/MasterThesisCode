#pragma once

#include <vector>
#include <array>
#include <complex>
#include <tuple>

template <class T, unsigned int H, unsigned int W = H>
class Matrix
{
public:
    static constexpr Matrix<T, H, H> Identity();

    Matrix();
    Matrix(std::initializer_list<std::initializer_list<T>> data);
    Matrix(std::vector<std::vector<T>> data);
    Matrix(std::vector<T> data);
    template <class ... F>
    Matrix(F... f);

    T& operator [] (unsigned int i);
    T& operator () (unsigned int i);
    T& operator () (unsigned int i, unsigned int j);

    operator T*() const;
    operator std::vector<T>() const;

    std::tuple<unsigned int, unsigned int> constexpr getSize();
    unsigned int constexpr getTotalSize();

    Matrix<T, H, W> operator - ();

    template <class F>
    Matrix<T, H, W> operator + (Matrix<F, H, W> m);
    template <class F>
    Matrix<T, H, W> operator - (Matrix<F, H, W> m);
    template <class F, unsigned int N>
    Matrix<T, H, N> operator * (Matrix<F, W, N> m);
    /*template <class F, typename = std::enable_if<W == 1>>
    T operator * (Matrix<F, H, 1> m)
    {
        T t = (*this)[0] * m[0];
        for (int i = 1; i < H; i++)
            t += (*this)[i] * m[i];
        return t;
    };*/
    template <class F>
    Matrix<T, H, W>& operator = (Matrix<F, H, W> m);
    template <class F>
    Matrix<T, H, W>& operator += (Matrix<F, H, W> m);
    template <class F>
    Matrix<T, H, W>& operator -= (Matrix<F, H, W> m);
    template <class F, typename = std::enable_if<W == H>>
    Matrix<T, H, W>& operator *= (Matrix<F, H, W> m)
    {
        *this = *this * m;
        return *this;
    };

    Matrix<T, H, W> operator + (T f);
    Matrix<T, H, W> operator - (T f);
    Matrix<T, H, W> operator * (T f);
    Matrix<T, H, W> operator / (T f);
    Matrix<T, H, W>& operator += (T f);
    Matrix<T, H, W>& operator -= (T f);
    Matrix<T, H, W>& operator *= (T f);
    Matrix<T, H, W>& operator /= (T f);

    template <unsigned int height, unsigned int width>
    Matrix<T, height, width> cut();
    Matrix<T, W, H> transpose();
    Matrix<T, W, H> adjoint();
protected:
    std::array<T, H * W> content;
};

template <class T, unsigned int H, unsigned int W>
Matrix<T, H, W> operator + (T f, Matrix<T, H, W> m);
template <class T, unsigned int H, unsigned int W>
Matrix<T, H, W> operator - (T f, Matrix<T, H, W> m);
template <class T, unsigned int H, unsigned int W>
Matrix<T, H, W> operator * (T f, Matrix<T, H, W> m);

template <class T, unsigned int N>
using Vector = Matrix<T, N, 1>;

template <class T, unsigned int H, unsigned int W>
Matrix<T, W, H> transpose(Matrix<T, H, W> m);
template <class T, unsigned int H, unsigned int W>
Matrix<T, W, H> adjoint(Matrix<T, H, W> m);


template<class T, unsigned int H, unsigned int W>
inline constexpr Matrix<T, H, H> Matrix<T, H, W>::Identity()
{
    Matrix<T, H, H> matrix;
    for (int i = 0; i < H; i++)
        matrix(i, i) = 1;
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>::Matrix()
{
    //content.fill(T(0));
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>::Matrix(std::initializer_list<std::initializer_list<T>> data)
{
    unsigned int i = 0, j = 0;
    for (auto row : data)
    {
        for (auto element : row)
        {
            content[i * W + j] = element;
            j++;
        }
        j = 0;
        i++;
    }
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>::Matrix(std::vector<std::vector<T>> data)
{
    if (data.size() != H)
        throw std::runtime_error("Data size doesn't match matrix size");
    for (unsigned int i = 0; i < H; i++)
    {
        if (data[i].size() != W)
            throw std::runtime_error("Data size doesn't match matrix size");
        for (unsigned int j = 0; j < W; j++)
        {
            content[i * W + j] = data[i][j];
        }
    }
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>::Matrix(std::vector<T> data)
{
    if (data.size() != H * W)
        throw std::runtime_error("Data size doesn't match matrix size");
    for (unsigned int i = 0; i < H * W; i++)
    {
        content[i] = data[i];
    }
}

template<class T, unsigned int H, unsigned int W>
template<class ...F>
inline Matrix<T, H, W>::Matrix(F ...f)
{
    T args[] = { static_cast<T>(std::forward<F>(f), f) ... };
    static_assert(H * W == sizeof...(F), "Argument count doesn't match vector size");
    for (unsigned int i = 0; i < sizeof...(F); i++)
        content[i] = args[i];
}

template<class T, unsigned int H, unsigned int W>
T& Matrix<T, H, W>::operator [] (unsigned int i)
{
    return content[i];
}

template<class T, unsigned int H, unsigned int W>
inline T& Matrix<T, H, W>::operator () (unsigned int i)
{
    return content[i];
}

template<class T, unsigned int H, unsigned int W>
inline T& Matrix<T, H, W>::operator () (unsigned int i, unsigned int j)
{
    return content[i * W + j];
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>::operator T*() const
{
    return &content[0];
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>::operator std::vector<T>() const
{
    return std::vector<T>(content.begin(), content.end());
}

template<class T, unsigned int H, unsigned int W>
inline constexpr std::tuple<unsigned int, unsigned int> Matrix<T, H, W>::getSize()
{
    return std::tuple<unsigned int, unsigned int>(H, W);
}

template<class T, unsigned int H, unsigned int W>
inline constexpr unsigned int Matrix<T, H, W>::getTotalSize()
{
    return H * W;
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W> Matrix<T, H, W>::operator - ()
{
    Matrix<T, H, W> matrix;
    for (int i = 0; i < H * W; i++)
        matrix[i] = -content[i];
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
template<class F>
inline Matrix<T, H, W> Matrix<T, H, W>::operator + (Matrix<F, H, W> m)
{
    Matrix<T, H, W> matrix;
    for (int i = 0; i < H * W; i++)
        matrix[i] = content[i] + m[i];
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
template<class F>
inline Matrix<T, H, W> Matrix<T, H, W>::operator - (Matrix<F, H, W> m)
{
    Matrix<T, H, W> matrix;
    for (int i = 0; i < H * W; i++)
        matrix[i] = content[i] - m[i];
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
template<class F, unsigned int N>
inline Matrix<T, H, N> Matrix<T, H, W>::operator * (Matrix<F, W, N> m)
{
    Matrix<T, H, N> matrix;
    for (unsigned int i = 0; i < H; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            matrix(i, j) = (*this)(i, 0) * m(0, j);
        }
    }
    for (unsigned int i = 0; i < H; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            for (unsigned int k = 1; k < W; k++)
            {
                matrix(i, j) += (*this)(i, k) * m(k, j);
            }
        }
    }
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
template<class F>
inline Matrix<T, H, W>& Matrix<T, H, W>::operator = (Matrix<F, H, W> m)
{
    for (int i = 0; i < H * W; i++)
        content[i] = m[i];
    return *this;
}

template<class T, unsigned int H, unsigned int W>
template<class F>
inline Matrix<T, H, W>& Matrix<T, H, W>::operator += (Matrix<F, H, W> m)
{
    for (int i = 0; i < H * W; i++)
        content[i] += m[i];
    return *this;
}

template<class T, unsigned int H, unsigned int W>
template<class F>
inline Matrix<T, H, W>& Matrix<T, H, W>::operator -= (Matrix<F, H, W> m)
{
    for (int i = 0; i < H * W; i++)
        content[i] -= m[i];
    return *this;
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W> Matrix<T, H, W>::operator + (T f)
{
    Matrix<T, H, W> matrix = *this;
    for (int i = 0; i < (std::min)(H, W); i++)
    {
        matrix(i, i) += f;
    }
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W> Matrix<T, H, W>::operator - (T f)
{
    Matrix<T, H, W> matrix = *this;
    for (int i = 0; i < (std::min)(H, W); i++)
    {
        matrix(i, i) -= f;
    }
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W> Matrix<T, H, W>::operator * (T f)
{
    Matrix<T, H, W> matrix;
    for (int i = 0; i < H * W; i++)
        matrix[i] = content[i] * f;
    return matrix;
}
template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W> Matrix<T, H, W>::operator / (T f)
{
    Matrix<T, H, W> matrix;
    for (int i = 0; i < H * W; i++)
        matrix[i] = content[i] / f;
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>& Matrix<T, H, W>::operator += (T f)
{
    *this = *this + f;
    return *this;
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>& Matrix<T, H, W>::operator -= (T f)
{
    *this = *this - f;
    return *this;
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>& Matrix<T, H, W>::operator *= (T f)
{
    for (int i = 0; i < H * W; i++)
        content[i] *= f;
    return *this;
}
template<class T, unsigned int H, unsigned int W>
inline Matrix<T, H, W>& Matrix<T, H, W>::operator /= (T f)
{
    for (int i = 0; i < H * W; i++)
        content[i] /= f;
    return *this;
}

template<class T, unsigned int H, unsigned int W>
template<unsigned int height, unsigned int width>
inline Matrix<T, height, width> Matrix<T, H, W>::cut()
{
    Matrix<T, height, width> matrix;
    unsigned int minHeight = (std::min)(height, H);
    unsigned int minWidth = (std::min)(width, W);
    for (int i = 0; i < minHeight; i++)
    {
        for (int j = 0; j < minWidth; j++)
        {
            matrix(i, j) = (*this)(i, j);
        }
    }
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, W, H> Matrix<T, H, W>::transpose()
{
    return transpose(*this);
}

template<class T, unsigned int H, unsigned int W>
inline Matrix<T, W, H> Matrix<T, H, W>::adjoint()
{
    return adjoint(*this);
}

template<class T, unsigned int H, unsigned int W>
Matrix<T, H, W> operator + (T f, Matrix<T, H, W> m)
{
    return m + f;
}

template<class T, unsigned int H, unsigned int W>
Matrix<T, H, W> operator - (T f, Matrix<T, H, W> m)
{
    return -m + f;
}

template<class T, unsigned int H, unsigned int W>
Matrix<T, H, W> operator * (T f, Matrix<T, H, W> m)
{
    return m * f;
}

template<class T, unsigned int H, unsigned int W>
Matrix<T, W, H> transpose(Matrix<T, H, W> m)
{
    Matrix<std::complex<T>, W, H> matrix;
    for (int i = 0; i < W; i++)
    {
        for (int j = 0; j < H; j++)
        {
            matrix(i, j) = m(j, i);
        }
    }
    return matrix;
}

template<class T, unsigned int H, unsigned int W>
Matrix<T, W, H> adjoint(Matrix<T, H, W> m)
{
    return transpose(m);
}

template<class T, unsigned int H, unsigned int W>
Matrix<std::complex<T>, W, H> adjoint(Matrix<std::complex<T>, H, W> m)
{
    Matrix<std::complex<T>, W, H> matrix;
    for (int i = 0; i < W; i++)
    {
        for (int j = 0; j < H; j++)
        {
            matrix(i, j) = std::conj(m(j, i));
        }
    }
    return matrix;
}
