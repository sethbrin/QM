#include <cstring>

#include "grid_structures.h"
#include "common.h"

using namespace database;

/**
 * DataStructure
 */
void DataStructure::set_param(const char* fftype)
{
  // set params
  DataStructureParam params_tmp;
  params_tmp.dim2[0] = 0;
  params_tmp.dim2[1] = 3.8210;
  m_params.insert(std::make_pair("CT-OW", params_tmp));

  params_tmp.dim2[1] = 2.2993;
  m_params.insert(std::make_pair("CT-HW", params_tmp));

  params_tmp.dim2[1] = 3.8279;
  m_params.insert(std::make_pair("HC-OW", params_tmp));

  params_tmp.dim2[1] = 1.4165;
  m_params.insert(std::make_pair("HC-HW", params_tmp));

  params_tmp.dim1 = 0;

  if (strcmp(fftype, "none") == 0)
  {
    m_params.insert(std::make_pair("qCT", params_tmp));
    m_params.insert(std::make_pair("qHC", params_tmp));
    m_params.insert(std::make_pair("qOw", params_tmp));
    m_params.insert(std::make_pair("qHW", params_tmp));
  }
  else
  {
    m_params.insert(std::make_pair("qCT", params_tmp));
    m_params.insert(std::make_pair("qHC", params_tmp));
    m_params.insert(std::make_pair("qOw", params_tmp));
    m_params.insert(std::make_pair("qHW", params_tmp));
  }

}

void DataStructure::set_symmetry()
{
  m_symface = {"xy"};
}

void DataStructure::set_num_of_atoms()
{
  m_n1 = 12;
  m_n2 = 3;
}

void DataStructure::set_R()
{
  m_R_NDX = {2.0,2.2,2.5,2.7,3.0,3.2,3.5,3.7,4.0,4.5,4.7,5.0,5.5,6.0,6.5,7.0,8.0};
  int size = m_R_NDX.size();

  for (int i=0; i<size-1; i++)
  {
    m_DR.push_back(m_R_NDX[i+1] - m_R_NDX[i]);
  }
  m_DR.push_back(0);
}

void DataStructure::set_phi()
{
  m_PHI_angles = {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0};
}

void DataStructure::set_theta()
{
  m_THETA_angles = {
    {0, {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0, 195.0, 210.0, 225.0, 240.0, 255.0, 270.0, 285.0, 300.0, 315.0, 330.0, 345.0}},
    {1, {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0, 195.0, 210.0, 225.0, 240.0, 255.0, 270.0, 285.0, 300.0, 315.0, 330.0, 345.0}},
    {2, {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0, 195.0, 210.0, 225.0, 240.0, 255.0, 270.0, 285.0, 300.0, 315.0, 330.0, 345.0}},
    {3, {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0}},
    {4, {0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0}},
    {5, {0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0}},
    {6, {0.0}}
  };
}

void DataStructure::set_nConf()
{
  m_nConf = 26;
  m_nNorm = 2;
}

void DataStructure::degree2radius()
{
  m_nGrid = 0;
  for (size_t i=0; i<m_PHI_angles.size(); i++)
  {
    m_PHI_angles[i] *= D2R;

    int THETA_anglesI = m_THETA_angles[i].size();
    m_NTheta.push_back(THETA_anglesI);
    for (int j=0; j<THETA_anglesI; j++)
    {
      m_THETA_angles[i][j] *= D2R;
    }
    m_nGrid += THETA_anglesI;
  }
}

void DataStructure::set_H_correction(bool flag)
{
  m_IsHmm = flag;
}

double DataStructure::lj_12_6(double r, double eps, double sig)
{
  double tmp = std::pow(sig/r, 6);

  return 4 * eps * (tmp * tmp - tmp);
}

double DataStructure::lj_9_6(double r, double eps, double sig)
{
  double tmp = std::pow(sig/r, 3);

  return eps * (2 * tmp * tmp * tmp - 3 * tmp * tmp);
}

double DataStructure::lj_b_14_7(double r, double eps, double sig)
{
  double rou = r / sig;
  return eps * std::pow(1.07 / (rou+0.07), 7) * (1.12/(std::pow(rou, 7)+0.12)-2.0);
}

double DataStructure::coulomb(double q0, double q1, double r)
{
  return q0 * q1 / r * 332.5;
}

double DataStructure::nocorr(double , double , double)
{
  return 0.0;
}

void DataStructure::initialize()
{
  set_param(m_fftype);
  set_symmetry();
  set_R();
  set_phi();
  set_theta();
  set_nConf();
  degree2radius();
  set_H_correction(false);
  set_num_of_atoms();

  std::map<const char*, Fun_lj> tmp = {
    {"none", &DataStructure::nocorr},
    {"9-6", &DataStructure::lj_9_6},
    {"12-6", &DataStructure::lj_12_6},
    {"14-7", &DataStructure::lj_b_14_7}
  };

  m_lj = tmp[m_fftype];

}



DataStructure::DataStructure(const char* fftype): m_fftype(fftype)
{
}

// Set the grid data of the structure_type
void DataStructure::set_grid_data(const char* structure_type)
{
  if (strcmp(structure_type, "wtr") == 0)
  {
    m_grid_data = {
         {0, { 0, 1, 13, 14}} ,
         {1, { 1, 2, 14, 15}} ,
         {2, { 2, 3, 15, 16}} ,
         {3, { 3, 4, 16, 17}} ,
         {4, { 4, 5, 17, 18}} ,
         {5, { 5, 6, 18, 19}} ,
         {6, { 6, 7, 19, 20}} ,
         {7, { 7, 8, 20, 21}} ,
         {8, { 8, 9, 21,  22}} ,
         {9, { 9, 10, 22, 23}} ,
        {10, {10, 11, 23, 24}} ,
        {11, {11, 12, 24, 25}} ,
        {13, {13, 14, 26, 27}} ,
        {14, {14, 15, 27, 28}} ,
        {15, {15, 16, 28, 29}} ,
        {16, {16, 17, 29, 30}} ,
        {17, {17, 18, 30, 31}} ,
        {18, {18, 19, 31, 32}} ,
        {19, {19, 20, 32, 33}} ,
        {20, {20, 21, 33, 34}} ,
        {21, {21, 22, 34, 35}} ,
        {22, {22, 23, 35, 36}} ,
        {23, {23, 24, 36, 37}} ,
        {24, {24, 25, 37, 38}} ,
        {26, {26, 27, 39, 40}} ,
        {27, {27, 28, 40, 41}} ,
        {28, {28, 29, 41, 41}} ,
        {29, {29, 30, 41, 42}} ,
        {30, {30, 31, 42, 43}} ,
        {31, {31, 32, 43, 44}} ,
        {32, {32, 33, 44, 44}} ,
        {33, {33, 34, 44, 45}} ,
        {34, {34, 35, 45, 46}} ,
        {35, {35, 36, 46, 47}} ,
        {36, {36, 37, 47, 47}} ,
        {37, {37, 38, 47, 48}} ,
        {39, {39, 40, 49, 50}} ,
        {40, {40, 41, 50, 50}} ,
        {41, {41, 42, 50, 51}} ,
        {42, {42, 43, 51, 52}} ,
        {43, {43, 44, 52, 52}} ,
        {44, {44, 45, 52, 53}} ,
        {45, {45, 46, 53, 54}} ,
        {46, {46, 47, 54, 54}} ,
        {47, {47, 48, 54, 55}} ,
        {49, {49, 50, 56, 57}} ,
        {50, {50, 51, 57, 57}} ,
        {51, {51, 52, 57, 58}} ,
        {52, {52, 53, 58, 59}} ,
        {53, {53, 54, 59, 59}} ,
        {54, {54, 55, 59, 60}} ,
        {56, {56, 57, 61, 61}} ,
        {57, {57, 58, 61, 61}} ,
        {58, {58, 59, 61, 61}} ,
        {59, {59, 60, 61, 61}}};

  }
}

std::pair<std::vector<double>, std::vector<double> > DataStructure::calt_vec1(const std::vector<double>& a0, const std::vector<double>& a1, const std::vector<double>& a2)
{
  assert(a0.size() == a1.size());
  assert(a0.size() == a2.size());
  return std::pair<std::vector<double>, std::vector<double> >(a0, {1, 0, 0});
}

std::pair<std::vector<double>, std::vector<double> > DataStructure::calt_vec2(const std::vector<double>& a0, const std::vector<double>& a1, const std::vector<double>& a2)
{
  assert(a0.size() == a1.size());
  assert(a0.size() == a2.size());
  return std::pair<std::vector<double>, std::vector<double> >(a1, {0, 1, 0});
}

std::vector<double> DataStructure::calt_dvec(const std::vector<double>& a0, const std::vector<double>& a1, const std::vector<double>& a2)
{
  assert(a0.size() == a1.size());
  assert(a0.size() == a2.size());
  double ratio_CA1 = 0.5618541311379575;
  double ratio_CA2 = 0.4381458688620425;

  int size = a0.size();
  std::vector<double> res(0, size);

  for (int i=0; i<size; i++)
  {
    res[i] = a0[i] * ratio_CA1 - a1[i] * ratio_CA2;
  }
  return res;
}


/**
 * PrpStructure
 */
void PrpStructure::set_param(const char* fftype)
{
  // set params
  DataStructureParam params_tmp;
  params_tmp.dim2[0] = 0;
  params_tmp.dim2[1] = 3.8210;
  m_params.insert(std::make_pair("CT-OW", params_tmp));

  params_tmp.dim2[1] = 2.2993;
  m_params.insert(std::make_pair("CT-HW", params_tmp));

  params_tmp.dim2[0] = 0.0100;
  params_tmp.dim2[1] = 2.5279;
  m_params.insert(std::make_pair("HC-OW", params_tmp));

  params_tmp.dim2[1] = 2.2165;
  m_params.insert(std::make_pair("HC-HW", params_tmp));

  params_tmp.dim1 = 0;

  if (strcmp(fftype,"none") == 0)
  {
    m_params.insert(std::make_pair("qCT", params_tmp));
    m_params.insert(std::make_pair("qHC", params_tmp));
    m_params.insert(std::make_pair("qOw", params_tmp));
    m_params.insert(std::make_pair("qHW", params_tmp));
  }
  else
  {
    m_params.insert(std::make_pair("qCT", params_tmp));
    m_params.insert(std::make_pair("qHC", params_tmp));
    params_tmp.dim1 = -0.6800;
    m_params.insert(std::make_pair("qOw", params_tmp));
    params_tmp.dim1 = 0.3400;
    m_params.insert(std::make_pair("qHW", params_tmp));
  }

}

void PrpStructure::set_theta()
{
  m_THETA_angles = {
    {0, {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0}},
    {1, {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0}},
    {2, {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0}},
    {3, {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0}},
    {4, {0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0}},
    {5, {0.0, 60.0, 120.0, 180.0}},
    {6, {0.0}}
  };
}

void PrpStructure::set_symmetry()
{
  m_symface = {"xy", "xz"};
}

void PrpStructure::set_num_of_atoms()
{
  m_n1 = 11;
  m_n2 = 3;
}

void PrpStructure::set_H_correction(bool flag)
{
  m_IsHmm = false;
  m_Hmm = {1,2,3,8,9,10};
}

std::pair<std::vector<double>, std::vector<double> > PrpStructure::calt_vec1(const std::vector<double>& a0, const std::vector<double>& a1, const std::vector<double>& a2)
{
  assert(a0.size() == a1.size());
  assert(a0.size() == a2.size());
  return std::pair<std::vector<double>, std::vector<double> >(a1, {1, 0, 0});
}

std::pair<std::vector<double>, std::vector<double> > PrpStructure::calt_vec2(const std::vector<double>& a0, const std::vector<double>& a1, const std::vector<double>& a2)
{
  assert(a0.size() == a1.size());
  assert(a0.size() == a2.size());

  int size = a0.size();
  std::vector<double> tmpa0(a0.size(), 0);
  std::vector<double> tmpa2(a0.size(), 0);

  for (int i=0; i<size; i++)
  {
    tmpa0[i] = a0[i] - a1[i];
    tmpa2[i] = a2[i] - a1[i];
  }

  return std::pair<std::vector<double>, std::vector<double> >(get_normal(tmpa0, tmpa2), {0,0,1});
}

std::vector<double> PrpStructure::calt_dvec(const std::vector<double>& a0, const std::vector<double>& a1, const std::vector<double>& a2)
{
  assert(a0.size() == a1.size());
  assert(a0.size() == a2.size());
  int size = a0.size();
  std::vector<double> res(0, size);

  for (int i=0; i<size; i++)
  {
    res[i] = -0.5*(0.5*(a0[i] + a2[i]) + a1[i]);
  }
  return res;
}


/**
 * WaterStructure
 */
void WaterStructure::set_num_of_atoms()
{
  m_n1 = 3;
  m_n2 = 3;
}

void WaterStructure::set_H_correction(bool flag)
{
  m_IsHmm = false;
  m_Hmm = {1,2,3,8,9,10};
}

void WaterStructure::set_theta()
{
  m_THETA_angles = {
    {0, {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0}},
    {1, {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0}},
    {2, {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0}},
    {3, {0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0}},
    {4, {0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0}},
    {5, {0.0, 45.0, 90.0, 135.0, 180.0}},
    {6, {0.0}}
  };
}

std::pair<std::vector<double>, std::vector<double> > WaterStructure::calt_vec1(const std::vector<double>& a0, const std::vector<double>& a1, const std::vector<double>& a2)
{
  assert(a0.size() == a1.size());
  assert(a0.size() == a2.size());
  int size = a0.size();

  std::vector<double> tmp(a0.size(), 0);
  for (int i=0; i<size; i++)
  {
    tmp[i] = (a1[i] + a2[i]) * 0.5;
  }
  return std::pair<std::vector<double>, std::vector<double> >(tmp, {1, 0, 0});
}

std::pair<std::vector<double>, std::vector<double> > WaterStructure::calt_vec2(const std::vector<double>& a0, const std::vector<double>& a1, const std::vector<double>& a2)
{
  assert(a0.size() == a1.size());
  assert(a0.size() == a2.size());

  return std::pair<std::vector<double>, std::vector<double> >(get_normal(a2, a1), {0,0,1});
}

std::vector<double> WaterStructure::calt_dvec(const std::vector<double>& a0, const std::vector<double>& a1, const std::vector<double>& a2)
{
  assert(a0.size() == a1.size());
  assert(a0.size() == a2.size());
  int size = a0.size();
  std::vector<double> res(size, 0);

  for (int i=0; i<size; i++)
  {
    res[i] = -a0[i];
  }
  return res;
}



// TODO this is an ugly design, later we should change it
// As if we set the fragtype first time, we can not change it
// it will last for the program end
DataStructure* database::getDataStructure(const char* fragtype)
{
  static DataStructure* pDS = NULL;
  if (pDS == NULL)
  {
    if (strcmp(fragtype,"wtr") == 0)
    {
      pDS = new WaterStructure(fragtype);
    } else
    {
      pDS = new DataStructure(fragtype);
    }
    pDS->set_grid_data(fragtype);
    pDS->initialize();
  }
  return pDS;
}


