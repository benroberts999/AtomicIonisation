#pragma once
#include "VectorOverloads.hpp"
#include <cassert>
#include <vector>

//! Basic class to store matrix
template <typename T = double> class Matrix {

protected:
  std::size_t m_rows;
  std::size_t m_cols;
  std::vector<T> m_data;

public:
  //! Initialise a blank matrix rows*cols, filled with 0
  Matrix(std::size_t rows = 0, std::size_t cols = 0)
      : m_rows(rows), m_cols(cols), m_data(rows * cols) {}

  //! Return rows [major index size]
  std::size_t rows() const { return m_rows; }

  //! Return columns [minor index size]
  std::size_t cols() const { return m_cols; }

  //! [] index access (with no range checking). [i][j]: ith row, jth col
  T *operator[](std::size_t i) { return &(m_data[i * m_cols]); }

  //! const [] index access (with no range checking). [i][j]: ith row, jth col
  const T *operator[](std::size_t i) const { return &(m_data[i * m_cols]); }

  //! () index access (with range checking) .at(i,j) returns ith row, jth col
  T &at(std::size_t row_i, std::size_t col_j) {
    assert(row_i < m_rows && col_j < m_cols);
    return m_data[row_i * m_cols + col_j];
  }

  //! const () index access (with range checking). at(i,j): ith row, jth col
  T at(std::size_t row_i, std::size_t col_j) const {
    assert(row_i < m_rows && col_j < m_cols);
    return m_data[row_i * m_cols + col_j];
  }

  //! () index access (with range checking). (i,j): ith row, jth col
  T &operator()(std::size_t i, std::size_t j) { return at(i, j); }

  //! const () index access (with range checking). (i,j): ith row, jth col
  T operator()(std::size_t i, std::size_t j) const { return at(i, j); }

  //============================================================================
  // Operator overloads: +,-, scalar */

  //! Provides += operator for two matrices (in-place adition)
  Matrix<T> &operator+=(const Matrix<T> &rhs) {
    assert(this->rows() == rhs.rows() && this->cols() == rhs.cols() &&
           "Matrices must have same dimensions for addition");
    using namespace VectorOverloads;
    this->m_data += rhs.m_data;
    return *this;
  }

  //! Provides -= operator for two matrices (in-place subtraction)
  Matrix<T> &operator-=(const Matrix<T> &rhs) {
    assert(rows() == rhs.rows() && cols() == rhs.cols() &&
           "Matrices must have same dimensions for subtraction");
    using namespace VectorOverloads;
    this->m_data -= rhs.m_data;
    return *this;
  }

  //! Provides in-place scalar multiplication ( M *= x )
  Matrix<T> &operator*=(const T x) {
    using namespace VectorOverloads;
    this->m_data *= x;
    return *this;
  }

  //============================================================================

  //! Provides matrix addition
  [[nodiscard]] friend Matrix<T> operator+(Matrix<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs += rhs);
  }

  //! Provides matrix subtraction
  [[nodiscard]] friend Matrix<T> operator-(Matrix<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs -= rhs);
  }

  //! Provides scalar multiplication ( x * M )
  [[nodiscard]] friend Matrix<T> operator*(const T x, Matrix<T> rhs) {
    return (rhs *= x);
  }

  //! Provides scalar multiplication ( M * x )
  [[nodiscard]] friend Matrix<T> operator*(Matrix<T> lhs, const T x) {
    return (lhs *= x);
  }
};