#pragma once
/**
 * \file common.h
 * \brief contain the helper functions
 */

#include <cstddef>
#include <cassert>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <functional>

using std::vector;

namespace common {

#define QMException(x)         \
  fprintf(stderr, "%s\n", #x); \
  exit(1);

vector<double> linspace(double min, double max, size_t n) {
  if (n == 1) return {min};
  double step = (max - min) / (n - 1);

  vector<double> res;
  double start = min;
  for (size_t idx = 0; idx < n; ++idx) {
    res.push_back(start);
    start += step;
  }
  return res;
}

vector<int> range(size_t n) {
  vector<int> res;
  for (size_t idx = 0; idx < n; ++idx) {
    res.push_back(idx);
  }
  return res;
}

////////////////////////// Vecop ///////////////////////
template <typename T>
vector<T> vec_op(const vector<T>& a, const vector<T>& b,
                 std::function<T(T, T)> op) {
  assert(a.size() == b.size());

  vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result),
                 op);
  return result;
}

template <typename T>
vector<T> vec_op(const vector<T>& a, T b,
                 std::function<T(T, T)> op) {
  vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), std::back_inserter(result),
                 std::bind(op, std::placeholders::_1, b));
  return result;
}

#define VECOP(op, opfun)                                          \
  template <typename T>                                           \
  vector<T> operator op(const vector<T>& a, const vector<T>& b) { \
    return vec_op<T>(a, b, std::opfun<T>());                      \
  }                                                               \
  template <typename T>                                           \
  vector<T> operator op(const vector<T>& a, T b) {                \
    return vec_op<T>(a, b, std::opfun<T>());                      \
  }                                                               \
  template <typename T>                                           \
  vector<T> operator op(T b, const vector<T>& a) {                \
    return vec_op<T>(a, b, std::opfun<T>());                      \
  }

VECOP(+, plus)
VECOP(-, minus)
VECOP(*, multiplies)
VECOP(/, divides)

}  // common
