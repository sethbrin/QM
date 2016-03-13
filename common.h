#ifndef _COMMON_H_
#define _COMMON_H_

#include <string>
#include <sstream>
#include <vector>

using std::vector;

class Coord
{
public:
  explicit Coord():
    m_x(0), m_y(0), m_z(0)
  {}
  explicit Coord(double x, double y, double z):
    m_x(x), m_y(y), m_z(z)
  {}
  double m_x;
  double m_y;
  double m_z;

  static double distance(Coord a, Coord b)
  {
    return (a.m_x-b.m_x)*(a.m_x - b.m_x)
      + (a.m_x-b.m_x)*(a.m_x - b.m_x)
      + (a.m_x-b.m_x)*(a.m_x - b.m_x);
  }
};

vector<double> translate(const vector<double>& vec, const vector<double>& dvec);
double dot(const vector<double>& v1, const vector<double>&v2);
vector<double> div(const vector<double>& vec, double factor);
vector<double> rotate(const vector<double>& vec, const vector<double> axis, double theta);
vector<double> get_normal(const vector<double>& vec1, const vector<double>& vec2);
double length(const vector<double>& vec);
double angle(const vector<double>& vec1, const vector<double>& vec2);
#endif
