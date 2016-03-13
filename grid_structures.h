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

    void initialize();

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
    const char* m_fftype;

    virtual void set_param(const char* fftype);
    virtual void set_symmetry();
    virtual void set_num_of_atoms();
    virtual void set_R();
    virtual void set_phi();
    virtual void set_theta();
    virtual void set_nConf();
    virtual void degree2radius();
    virtual void set_H_correction(bool flag);

    virtual std::pair<std::vector<double>, std::vector<double> > calt_vec1(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
    virtual std::pair<std::vector<double>, std::vector<double> > calt_vec2(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
    virtual std::vector<double> calt_dvec(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);

    double lj_12_6(double r, double eps, double sig);
    double lj_9_6(double r, double eps, double sig);
    double lj_b_14_7(double r, double eps, double sig);
    double nocorr(double r, double eps, double sig);
  private:
  };

  class PrpStructure: public DataStructure
  {
  public:
    explicit PrpStructure(const char* fftype):DataStructure(fftype){}
    void set_param(const char* fftype);
    void set_theta();
    void set_symmetry();
    void set_num_of_atoms();
    void set_H_correction(bool flag);

    std::pair<std::vector<double>, std::vector<double> > calt_vec1(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
    std::pair<std::vector<double>, std::vector<double> > calt_vec2(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
    std::vector<double> calt_dvec(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
  };

  class WaterStructure : public PrpStructure
  {
  public:
    explicit WaterStructure(const char* fftype):PrpStructure(fftype){}
    void set_num_of_atoms();
    void set_theta();
    void set_H_correction(bool flag);

    std::pair<std::vector<double>, std::vector<double> > calt_vec1(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
    std::pair<std::vector<double>, std::vector<double> > calt_vec2(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);
    std::vector<double> calt_dvec(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c);

    //virtual void set_num_of_atoms();
    static const std::map<int, std::vector<int> > grid_data;
  private:
  };


  DataStructure* getDataStructure(const char* fragtype);
}
#endif
