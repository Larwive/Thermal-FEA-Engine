#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <vector>

namespace thermal_fea::linalg {

class Vector {
public:
  Vector(std::size_t size) : data_(size, .0) {}

  std::size_t size() const;

  double &operator()(std::size_t i);

  const double &operator()(std::size_t i) const;

  double dot(const Vector &other) const;

  double norm() const;

  friend std::ostream &operator<<(std::ostream &os, const Vector &v);

private:
  std::vector<double> data_;
};

Vector operator+(const Vector &a, const Vector &b);

Vector operator-(const Vector &a, const Vector &b);

Vector operator*(double alpha, const Vector &v);

Vector &operator+=(Vector &a, const Vector &b);

} // namespace thermal_fea::linalg
