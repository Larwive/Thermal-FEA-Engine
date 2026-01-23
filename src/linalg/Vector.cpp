#include <cassert>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <thermal_fea/linalg/Vector.hpp>

namespace thermal_fea::linalg {

std::size_t Vector::size() const { return data_.size(); }

double &Vector::operator()(std::size_t i) { return data_[i]; }

const double &Vector::operator()(std::size_t i) const { return data_[i]; }

double Vector::dot(const Vector &other) const {
  assert(size() == other.size());
  double s = 0.0;
  for (std::size_t i = 0; i < size(); ++i)
    s += data_[i] * other.data_[i];
  return s;
}

double Vector::norm() const { return std::sqrt(dot(*this)); }

std::ostream &operator<<(std::ostream &os, const Vector &v) {
  for (std::size_t i = 0; i < v.size(); ++i) {
    os << v(i);
    if (i + 1 < v.size())
      os << " ";
  }
  return os;
}

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
