#include <vector>
#include <map>
#include <array>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "common.h"

using std::vector;
using std::pair;
using std::array;

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


vector<double> rotate(const vector<double>& vec, const vector<double>& axis, double theta)
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

  vector<double> res(rotmat.size(), 0);
  for(unsigned int i=0; i< rotmat.size();i++)
  {
    res[i] = dot(rotmat[i], vec);
  }
  return res;


}

double distance(const vector<double>& vec1, const vector<double>& vec2)
{
  assert(vec1.size() == 3);
  assert(vec2.size() == 3);
  double res = 0;
  res += vec1[1] * vec2[2] - vec2[1] * vec1[2];
  res += vec1[2] * vec2[0] - vec2[2] * vec1[0];
  res += vec1[0] * vec2[1] - vec2[1] * vec1[0];
  return res;

}

vector<double> get_normal(const vector<double>& vec1, const vector<double>& vec2)
{
  assert(vec1.size() == 3);
  assert(vec2.size() == 3);
  vector<double> res(3,0.0);
  res[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
  res[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
  res[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
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

vector<double> subtraction(const vector<double>& a, const vector<double>& b)
{
  assert(a.size() == b.size());
  vector<double> res(a.size(), 0);

  for (size_t i=0; i<a.size(); i++)
  {
    res[i] = a[i] - b[i];
  }

  return res;
}


vector<double> get_bisect_unit(const vector<double>& a, const vector<double>& b)
{
  assert(a.size() == 3);
  assert(b.size() == 3);
  vector<double> res(3, 0);
  for (size_t i=0; i<3; i++)
  {
    res[i] = (a[i] + b[i])/2;
  }

  get_unit(res);
  return res;
}

void get_unit(vector<double>& vec)
{

  double sum = 0;
  for (double& item: vec)
  {
    sum += item * item;
  }
  sum = ::sqrt(sum);
  for (double& item: vec)
  {
    item /= sum;
  }
}

vector<double> get_normal_unit(const vector<double>& a, const vector<double>& b)
{
  assert(a.size() == 3);
  assert(b.size() == 3);

  vector<double> res = {a[1]*b[2] - a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};

  get_unit(res);

  return res;
}

std::string static quadrant_ndx(const vector<double> vec)
{
  std::array<std::string, 3> ndx = {"1", "1", "1"};
  for (size_t i=0; i<3; i++)
  {
    if (vec[i] < 0)
    {
      ndx[i] = "0";
    }
  }

  return ndx[0] + ndx[1] + ndx[2];
}


double linear1(double dxx0, double dx1x0, double y0, double y1)
{
  return y0 + (y1-y0)*dxx0/dx1x0;
}

/**
 *
 * points = [(ang1,ener1),(ang2,ener2),(ang3,ener3)]
 * x: the query angle
 * return the energy at point x
 */
double lagrange_interp(const std::vector<std::pair<double, double>>& points, double x) {
    double total = 0;
    int n = points.size();

    auto g = [&points, x,n](int i) {
        double tot_mul = 1;
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            tot_mul *= (x- points[j].first) /(points[i].first - points[j].first);
        }
        return tot_mul;
    };
    for (int i = 0; i < n; i++){
        total += points[i].second * g(i);
    }
    return total;
}

double bilinear(double y0, double y1, double y2, double y3, double angx, double angy)
{
  /**
   *Bilinear:
   *
   *   y3<----y2
   *          ^
   *      ya  | ---
   *          |  |angy
   *   y0---->y1---
   *   |<->|
   *    angx
   */
  double pi4 = M_PI / 4;
  double y01 = linear1(angx, pi4, y0, y1);
  double y32 = linear1(angx, pi4, y3, y2);

  return linear1(angy, pi4, y01, y32);
}

double bilinear_gen(double y0, double y1, double y2, double y3, double angx1, double angx2, double angy, int label) {
    double y01 = linear1(angx1, 1, y0, y1);
    double y23;
    if (label == 1) {
        y23 = linear1(angx2, 1, y2, y3);
    } else {
        y23 = linear1(angx2, 1, y3, y2);
    }
    return linear1(angy, 1, y01, y23);
}

std::map<std::string, std::map<int, std::vector<int>>> Mol2BisectDataNdx ={
              {"111", {{0,{ 4,20,19, 3}}, {1,{25,18,19,20}}, {2,{ 2, 3,19,18}}}},
              {"100", {{0,{ 4,12,13, 5}}, {1,{24,14,13,12}}, {2,{ 6, 5,13,14}}}},
              {"010", {{0,{ 0,16,23, 7}}, {1,{25,22,23,16}}, {2,{ 6, 7,23,22}}}},
              {"001", {{0,{ 0, 8, 9, 1}}, {1,{24,10, 9, 8}}, {2,{ 2, 1, 9,10}}}},
              {"110", {{0,{ 4, 5,21,20}}, {1,{25,20,21,22}}, {2,{ 6,22,21, 5}}}},
              {"101", {{0,{ 4, 3,11,12}}, {1,{24,12,11,10}}, {2,{ 2,10,11, 3}}}},
              {"011", {{0,{ 0, 1,17,16}}, {1,{25,16,17,18}}, {2,{ 2,18,17, 1}}}},
              {"000", {{0,{ 0, 7,15, 8}}, {1,{24, 8,15,14}}, {2,{ 6,14,15, 7}}}}
};
// for nh4,ch4 probe
//Mol2BisectDataNdx = { '111': {0:[ 0, 8,15, 7], 1:[24,14,15, 8], 2:[ 6, 7,15,14]},
//              '100': {0:[ 0,16,17, 1], 1:[16,25,18,17], 2:[ 2, 1,17,18]},
//              '010': {0:[ 4,12,11, 3], 1:[24,10,11,12], 2:[ 2, 3,11,10]},
//              '001': {0:[ 4,20,21, 5], 1:[25,22,21,20], 2:[ 6, 5,21,22]},
//              '110': {0:[ 0, 1, 9, 8], 1:[24, 8, 9,10], 2:[ 2,10, 9, 1]},
//              '101': {0:[ 0, 7,23,16], 1:[25,16,23,22], 2:[ 6,22,23, 7]},
//              '011': {0:[ 4, 5,13,12], 1:[24,12,13,14], 2:[ 6,14,13, 5]},
//              '000': {0:[ 4, 3,19,20], 1:[25,20,19,18], 2:[ 2,18,19, 3]}}

std::map<std::string, std::vector<std::vector<double>>> AxisInQuadrant = {
              {"111", {{ 1, 0, 0},{ 0, 1, 0},{ 0, 0, 1}}},
              {"100", {{ 1, 0, 0},{ 0,-1, 0},{ 0, 0,-1}}},
              {"010", {{-1, 0, 0},{ 0, 1, 0},{ 0, 0,-1}}},
              {"001", {{-1, 0, 0},{ 0,-1, 0},{ 0, 0, 1}}},
              {"110", {{ 1, 0, 0},{ 0, 1, 0},{ 0, 0,-1}}},
              {"101", {{ 1, 0, 0},{ 0,-1, 0},{ 0, 0, 1}}},
              {"011", {{-1, 0, 0},{ 0, 1, 0},{ 0, 0, 1}}},
              {"000", {{-1, 0, 0},{ 0,-1, 0},{ 0, 0,-1}}}
};


vector<int> weights_in_subsection(const vector<double>& bisvec, double& wghx, double& wghy, double cutoff)
{
  std::string quad_ndx = quadrant_ndx(bisvec);

  std::vector<double> angcos = {
    dot(bisvec, AxisInQuadrant[quad_ndx][0]),
    dot(bisvec, AxisInQuadrant[quad_ndx][1]),
    dot(bisvec, AxisInQuadrant[quad_ndx][2])
  };

  double maxcos = angcos[0];
  int subndx = 0;

  //The angle with max cos is the smallest
  for (int i=1; i<=2; i++)
  {
    if (angcos[i] > maxcos)
    {
      maxcos = angcos[i];
      subndx = i;
    }
  }

  /*
   **
   * Bilinear interpolation is used
   *     for the interpolation in a sub-section of a quadrant:
   * Three SUBSECTIONS: below is a quadrant with three sub-sections:
   *      z
   *     / \
   *    /   \
   *   c     b
   *  /  \o/  \
   * /    |    \
   *x-----a-----y
   *
   *subsection 0: x-a-o-c (0-1-4-3)
   *           1: y-b-o-a
   *           2: z-c-o-b
   * indices of axis: 0, 1, 2, 3, 4 (3,4 for cyclic call)
  */
  std::map<std::string, std::map<int, std::vector<double>>> ax = {
        {"111",{{0,{0,1,2}}, {1,{1,2,0}}, {2,{2,0,1}}}},
        {"100",{{0,{0,1,2}}, {1,{1,2,0}}, {2,{2,0,1}}}},
        {"010",{{0,{0,1,2}}, {1,{1,2,0}}, {2,{2,0,1}}}},
        {"001",{{0,{0,1,2}}, {1,{1,2,0}}, {2,{2,0,1}}}},
        {"110",{{0,{0,2,1}}, {1,{1,0,2}}, {2,{2,1,0}}}},
        {"101",{{0,{0,2,1}}, {1,{1,0,2}}, {2,{2,1,0}}}},
        {"011",{{0,{0,2,1}}, {1,{1,0,2}}, {2,{2,1,0}}}},
        {"000",{{0,{0,2,1}}, {1,{1,0,2}}, {2,{2,1,0}}}}
  };

  // angx is calculated by arctan(dy/dx)
  double tempx = bisvec[ax[quad_ndx][subndx][0]];
  double tempy = bisvec[ax[quad_ndx][subndx][1]];
  double angx = std::atan(std::fabs(tempy / tempx));

  // angy is calculated by arctan(dy/dx):
  double tempz = bisvec[ax[quad_ndx][subndx][2]];
  double angy = std::asin(std::fabs(tempz));
  // Bilinear:
  //
  //  y3<----y2
  //         ^
  //     ya  | ---
  //         |  |angy
  //  y0---->y1---
  //  |<->|
  //   angx
  //
  // np.pi/4.0:
  //pi4 = 0.78539816339744817
  //
  //y01 = linear1( angx, pi4, y0, y1 )
  //y32 = linear1( angx, pi4, y3, y2 )
  //ya  = linear1( angy, pi4, y01, y32 )

  wghx = angx;
  wghy = angy;
  return {
    Mol2BisectDataNdx[quad_ndx][subndx][0],
    Mol2BisectDataNdx[quad_ndx][subndx][1],
    Mol2BisectDataNdx[quad_ndx][subndx][2],
    Mol2BisectDataNdx[quad_ndx][subndx][3]
  };

}

/**
 * Return the weights of the proper two configurations at one
 * bisector direction.
 * Linear interpolation is used
 */

void weights_for_normal_general(const vector<double>& normal_vec, const vector<vector<double>>& config_vecs, double& w1, double& w2, int& ndx1, int& ndx2, double cutoff) {
    int nvec = config_vecs.size();
    vector<double> vec1 = config_vecs[0];
    vector<double> vec2 = config_vecs[1];

    vector<double> vx = vec1;
    vector<double> vz = get_normal_unit(vec1, vec2);
    vector<double> vy = get_normal(vz, vx);

    double da = M_PI / nvec;

    double px = dot(normal_vec, vx);
    double py = dot(normal_vec, vy);

    double rrot = 0;
    if (std::abs(px) < cutoff) {
        if (py>0) {
            w2 = 1.0;
        } else {
            w2 = -1.0;
        }
        w1 = 0;
        ndx2 = static_cast<int>((M_PI/2)/(M_PI/nvec));
        ndx1 = ndx2 - 1;
    } else {
        double ang = atan(py/px);
        if (px < 0) {
            ang += M_PI;
        } else if (py < 0) {
            ang += 2 * M_PI;
        }
        rrot = ang - static_cast<int>(ang/M_PI) * M_PI;
    }

    ndx1 = static_cast<int>(rrot / (da));
    double ang2 = rrot - static_cast<int>(rrot/da) * da;
    w2 = ang2 / da;
    w1 = 1 - w2;
    if (ndx1 == nvec - 1) ndx2 = 0;
    else ndx2 = ndx1 + 1;
}

pair<double, array<int, 3>> get_neighors_for_normal(const vector<double>& normal_vec, const vector<vector<double>>& config_vecs, double cutoff) {
    int nvec = config_vecs.size();
    vector<double> vec1 = config_vecs[0];
    vector<double> vec2 = config_vecs[1];

    vector<double> vx = vec1;
    vector<double> vz = get_normal_unit(vec1, vec2);
    vector<double> vy = get_normal(vz, vx);

    double da = M_PI / nvec;

    double px = dot(normal_vec, vx);
    double py = dot(normal_vec, vy);

    double ang0 = 0;
    int ndx1;
    int ndx2;
    int ndx3;
    if (std::abs(px) < cutoff) {
        ang0 = M_PI / 2;
        ndx2 =  static_cast<int>(ang0 / da);
        ndx1 = ndx2 - 1;

        if (ndx2 == nvec - 1) {
            ndx3 = nvec;
            if (ndx1 == 0) ERROR("ndx1 == 0");
        } else {
            ndx3 = ndx2 + 1;
        }
    } else {
        double ang = atan(py/px);
        if (px < 0) {
            ang += M_PI;
        } else if (py < 0) {
            ang += 2 * M_PI;
        }
        double rrot = ang - static_cast<int>(ang/M_PI) * M_PI;
        ang0 = rrot;

        ndx1 = static_cast<int>(rrot / da);
        if (ndx1 == nvec-1) {
            ndx2 = nvec;
            ndx3 = ndx1 - 1;
            if (ndx3 == 0) ERROR("ndx3 == 0");
        } else {
            ndx2 = ndx1 + 1;
            double tmp1 = (ndx2 + 1) * da - ang0;
            double tmp2 = ang0 - (ndx1 - 1) * da;
            if (std::abs(tmp1) < std::abs(tmp2)) {
                if (ndx2 == nvec-1) {
                    ndx3 = nvec;
                } else {
                    ndx3 = ndx2 + 1;
                }
            } else {
                ndx3 = ndx1 - 1;
            }
        }
    }
    return {ang0, {ndx1, ndx2, ndx3}};
}
