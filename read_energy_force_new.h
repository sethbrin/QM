/*
 * read energe force data from database
 *
 * Author: Ping Zeng(zengping10212044@gmail.com)
 */
#ifndef _READ_ENERGY_FORCE_NEW_H_
#define _READ_ENERGY_FORCE_NEW_H_

#include <cmath>
#include <map>
#include <vector>
#include <cassert>

#include "grid_structures.h"

namespace database
{
  class Config
  {
  };

  class DataBase
  {
  };

  class EnergeForceDatabase
  {
  public:
    explicit EnergeForceDatabase(const char* filename, const char* fragtype);
    ~EnergeForceDatabase();
    Config& at_all(int dim1, int dim2, int dim3, int dim4);
    Config& at_wtr(int dim1, int dim2);
    void read_file();

  private:

    Config* m_all_config;
    Config* m_wtr_config;
    std::string m_filename;
    std::string m_fragtype;
    DataStructure* m_pDS;

  };
}

#endif
