#pragma once

#include <cmath>
#include <cstddef>
#include <unordered_map>
// #include <vector>

#include "thermal_fea/linalg/Vector.hpp"

namespace thermal_fea::linalg {

/*class Matrix {
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
  };*/

class Matrix {
public:
  Matrix(std::size_t rows, std::size_t cols) : rows_(rows), cols_(cols) {}

  class Entry {
  public:
    Entry(Matrix &A, std::size_t i, std::size_t j) : A_(A), i_(i), j_(j) {}

    Entry &operator=(double v) {
      A_.set(i_, j_, v);
      return *this;
    }

    Entry &operator+=(double v) {
      A_.set(i_, j_, A_.get(i_, j_) + v);
      return *this;
    }

    Entry &operator-=(double v) {
      A_.set(i_, j_, A_.get(i_, j_) - v);
      return *this;
    }

    operator double() const { return A_.get(i_, j_); }

  private:
    Matrix &A_;
    std::size_t i_, j_;
  };

  std::size_t rows() const;
  std::size_t cols() const;

  Entry operator()(std::size_t i, std::size_t j) { return Entry(*this, i, j); }

  double operator()(std::size_t i, std::size_t j) const;

private:
  static constexpr double eps = 1e-12;

  void set(std::size_t i, std::size_t j, double v);
  double get(std::size_t i, std::size_t j) const;
  void insert_or_update(std::size_t i, std::size_t j, double v);
  void erase(std::size_t i, std::size_t j);

  std::size_t rows_;
  std::size_t cols_;
  std::unordered_map<std::size_t, double> data_;
};

Vector operator*(const Matrix &A, const Vector &x);

Matrix operator*(const double a, const Matrix &A);

} // namespace thermal_fea::linalg
