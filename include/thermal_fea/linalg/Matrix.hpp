#pragma once

#include <cstddef>
#include <vector>

#include <thermal_fea/linalg/Vector.hpp>

namespace thermal_fea::linalg {

class Matrix {
public:
  Matrix(std::size_t rows, std::size_t cols)
      : rows_(rows), cols_(cols), data_(rows * cols, .0) {}

  size_t rows() const;

  size_t cols() const;

  double &operator()(std::size_t i, std::size_t j);

  const double &operator()(std::size_t i, std::size_t j) const;

private:
  std::size_t rows_;
  std::size_t cols_;
  std::vector<double> data_;
};

Vector operator*(const Matrix &A, const Vector &x);

Matrix operator*(const double a, const Matrix &A);

} // namespace thermal_fea::linalg
