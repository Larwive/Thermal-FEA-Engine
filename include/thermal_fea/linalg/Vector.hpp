#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

namespace thermal_fea::linalg {

// TODO: Move implementation to Vector.cpp
class Vector {
public:
  Vector(std::size_t size) : data_(size, .0) {}

  std::size_t size() const { return data_.size(); }

  double &operator()(std::size_t i) { return data_[i]; }

  const double &operator()(std::size_t i) const { return data_[i]; }

  double dot(const Vector &other) const {
    assert(size() == other.size());
    double s = 0.0;
    for (std::size_t i = 0; i < size(); ++i)
      s += data_[i] * other.data_[i];
    return s;
  }

  double norm() const { return std::sqrt(dot(*this)); }

private:
  std::vector<double> data_;
};

Vector operator+(const Vector &a, const Vector &b) {
  assert(a.size() == b.size());
  Vector r(a.size());
  for (std::size_t i = 0; i < a.size(); ++i)
    r(i) = a(i) + b(i);
  return r;
}

Vector operator-(const Vector &a, const Vector &b) {
  assert(a.size() == b.size());
  Vector r(a.size());
  for (std::size_t i = 0; i < a.size(); ++i)
    r(i) = a(i) - b(i);
  return r;
}

Vector operator*(double alpha, const Vector &v) {
  Vector r(v.size());
  for (std::size_t i = 0; i < v.size(); ++i)
    r(i) = alpha * v(i);
  return r;
}

Vector &operator+=(Vector &a, const Vector &b) {
  assert(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i)
    a(i) += b(i);
  return a;
}

} // namespace thermal_fea::linalg
