/*
 * define some constant of QM calc
 *
 * Author: Ping Zeng(zengping10212044@gmail.com)
 */
#ifndef _GRID_STRUCTURES_H_
#define _GRID_STRUCTURES_H_

#include <cmath>
#include <map>
#include <vector>
#include <cassert>

namespace database
{

  typedef union dataStructureParam
  {
    double dim1;
    double dim2[2];
  } DataStructureParam;

  const double D2R = 3.14159265358979/180.0;

  class DataStructure;
  typedef double (DataStructure::*Fun_lj)(double, double, double);

  class WaterStructure;

  class DataStructure
  {
  public:
    explicit DataStructure(const char* fftype);


    double coulomb(double q0, double q1, double r);

    std::map<const char*, DataStructureParam> m_params;
    std::vector<std::string> m_symface;
    int m_n1;
    int m_n2;
    std::vector<double> m_R_NDX;
    std::vector<double> m_DR;
    std::vector<double> m_PHI_angles;
    std::map<int, std::vector<double> > m_THETA_angles;
    int m_nConf;
    int m_nNorm;
    std::vector<int> m_NTheta;
    int m_nGrid;
    bool m_IsHmm;
    std::vector<int> m_Hmm;
    Fun_lj m_lj;

  private:
    void set_param(const char* fftype);
    void set_symmetry();
    virtual void set_num_of_atoms();
    void set_R();
    void set_phi();
    void set_theta();
    void set_nConf();
    void degree2radius();
    void set_H_correction(bool flag);
    double lj_12_6(double r, double eps, double sig);
    double lj_9_6(double r, double eps, double sig);
    double lj_b_14_7(double r, double eps, double sig);
    double nocorr(double r, double eps, double sig);
  };


  class WaterStructure : public DataStructure
  {
  public:
    explicit WaterStructure(const char* fftype);

    //virtual void set_num_of_atoms();
    static const std::map<int, std::vector<int> > grid_data;
  private:
  };


  DataStructure* getDataStructure(const char* fragtype);
}
#endif
