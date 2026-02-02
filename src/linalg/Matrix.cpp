#include <cstddef>
// #include <vector>

#include <omp.h>

#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/linalg/Vector.hpp"

namespace thermal_fea::linalg {

size_t Matrix::rows() const { return rows_; }

size_t Matrix::cols() const { return cols_; }

/*double &Matrix::operator()(std::size_t i, std::size_t j) {
  return data_[i * cols_ + j];
}*/

/*const double &Matrix::operator()(std::size_t i, std::size_t j) const {
  return data_[i * cols_ + j];
  }*/
static std::size_t key(std::size_t i, std::size_t j, std::size_t cols) {
  return i * cols + j;
}

double Matrix::get(std::size_t i, std::size_t j) const {
  auto it = data_.find(key(i, j, cols_));
  return it == data_.end() ? 0.0 : it->second;
}

void Matrix::set(std::size_t i, std::size_t j, double v) {
  if (std::abs(v) < eps)
    erase(i, j);
  else
    insert_or_update(i, j, v);
}

void Matrix::insert_or_update(std::size_t i, std::size_t j, double v) {
  data_[key(i, j, cols_)] = v;
}

void Matrix::erase(std::size_t i, std::size_t j) {
  data_.erase(key(i, j, cols_));
}

double Matrix::operator()(std::size_t i, std::size_t j) const {
  return get(i, j);
}

Vector operator*(const Matrix &A, const Vector &x) {
  assert(A.cols() == x.size());

  Vector y(A.rows());
#pragma omp parallel for
  for (std::size_t i = 0; i < A.rows(); ++i) {
    double sum = 0.0;
#pragma omp simd reduction(+:sum)
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
