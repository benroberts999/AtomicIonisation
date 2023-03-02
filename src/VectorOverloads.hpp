#pragma once
#include <vector>

namespace VectorOverloads {

//! In-place addition of two std::vectors. They don't need to be the same length
template <typename T>
std::vector<T> &operator+=(std::vector<T> &a, const std::vector<T> &b) {
  const auto size = std::max(a.size(), b.size());
  a.resize(size); // T must be default constructable
  for (auto i = 0ul; i < b.size(); ++i) {
    a[i] += b[i];
  }
  return a;
}

//! Addition of two std::vectors
template <typename T>
std::vector<T> operator+(std::vector<T> a, const std::vector<T> &b) {
  return a += b;
}

//! In-place subtraction of two std::vectors
template <typename T>
std::vector<T> &operator-=(std::vector<T> &a, const std::vector<T> &b) {
  const auto size = std::max(a.size(), b.size());
  a.resize(size); // T must be default constructable
  for (auto i = 0ul; i < b.size(); ++i) {
    a[i] -= b[i];
  }
  return a;
}

//! Subtraction of two std::vectors
template <typename T>
std::vector<T> operator-(std::vector<T> a, const std::vector<T> &b) {
  return a -= b;
}

//! In-place scalar multiplication of std::vector
template <typename T, typename U>
std::vector<T> &operator*=(std::vector<T> &v, U x) {
  for (auto &v_i : v) {
    v_i *= x;
  }
  return v;
}

//! Scalar multiplication of std::vector, ( v * x )
template <typename T, typename U>
std::vector<T> operator*(std::vector<T> v, U x) {
  return v *= x;
}

//! Scalar multiplication of std::vector, ( x * v )
template <typename T, typename U>
std::vector<T> operator*(U x, std::vector<T> v) {
  return v *= x;
}

} // namespace VectorOverloads