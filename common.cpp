#include <vector>
#include <cmath>
#include <cassert>

#include "common.h"

using std::vector;

vector<double> translate(const vector<double>& vec, const vector<double>& dvec)
{
  vector<double> res;
  assert(vec.size()>=3);
  assert(dvec.size()>=3);

  for(int i=0; i<3;i++)
    res.push_back(vec[i] + dvec[i]);
  return res;
}

//rotate 

// 向量的点乘函数
double dot(const vector<double>& v1, const vector<double>&v2)
{
  assert(v1.size() == v2.size());
  double res = 0.0;
  for(unsigned int i=0; i<v1.size();i++)
  {
    res += v1[i] * v2[i];
  }
  return res;

}

// 向量的除以标量
vector<double> div(const vector<double>& vec, double factor)
{
  assert(factor != 0);
  vector<double> res(vec);
  for(unsigned int i=0; i< vec.size();i++)
    res[i] = res[i]/factor;

  return res;
}


vector<double> rotate(const vector<double>& vec, const vector<double> axis, double theta)
{
  double temp = sqrt(dot(axis,axis));
  vector<double> newAxis = div(axis,temp);

  double a = cos(0.5 * theta);

  double s = sin(0.5 * theta);
  double b = newAxis[0] * s;
  double c = newAxis[1] * s;
  double d = newAxis[2] * s;

  vector<vector<double>> rotmat = {{a*a+b*b-c*c-d*d,2*(b*c-a*d),2*(b*d+a*c)},{2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)},
    {2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c}};

  vector<double> res;
  for(unsigned int i=0; i< rotmat.size();i++)
  {
    res.push_back(dot(rotmat[i],vec));
  }
  return res;


}

vector<double> get_normal(const vector<double>& vec1, const vector<double>& vec2)
{
  assert(vec1.size() == 3);
  assert(vec2.size() == 3);
  vector<double> res(3,0.0);
  res[0] = vec1[1] * vec2[2] - vec2[1] * vec1[2];
  res[1] = vec1[2] * vec2[0] - vec2[2] * vec1[0];
  res[2] = vec1[0] * vec2[1] - vec2[1] * vec1[0];
  return res;

}

double length(const vector<double>& vec)
{
  return sqrt(dot(vec,vec));
}

double angle(const vector<double>& vec1, const vector<double>& vec2)
{
  double dotprod = dot(vec1,vec2);
  return acos(dotprod/((length(vec1) * length(vec2))));
}
