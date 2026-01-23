#include <cstddef>
#include <vector>

#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/linalg/Vector.hpp"

namespace thermal_fea::linalg {

size_t Matrix::rows() const { return rows_; }

size_t Matrix::cols() const { return cols_; }

double &Matrix::operator()(std::size_t i, std::size_t j) {
  return data_[i * cols_ + j];
}

const double &Matrix::operator()(std::size_t i, std::size_t j) const {
  return data_[i * cols_ + j];
}

Vector operator*(const Matrix &A, const Vector &x) {
  assert(A.cols() == x.size());

  Vector y(A.rows());
  for (std::size_t i = 0; i < A.rows(); ++i) {
    double sum = 0.0;
    for (std::size_t j = 0; j < A.cols(); ++j) {
      sum += A(i, j) * x(j);
    }
    y(i) = sum;
  }
  return y;
}

Matrix operator*(const double a, const Matrix &A) {

  Matrix B(A.rows(), A.cols());

  for (std::size_t i = 0; i < A.rows(); ++i) {
    for (std::size_t j = 0; j < A.cols(); ++j) {
      B(i, j) = a * A(i, j);
    }
  }
  return B;
}

} // namespace thermal_fea::linalg
