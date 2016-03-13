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
#include "gzstream/gzstream.h"
#include "common.h"

namespace database
{
  class Atom
  {
  public:
    explicit Atom(std::string& line);

    std::string m_name;
    Coord m_coord;
  };


  class Config
  {
  public:
    explicit Config();
    explicit Config(int n1, int n2);

    void setEnergy(double energe);
    void setName(std::string name);
    int m_n1;
    int m_n2;
    int m_energy;
    std::string m_name;
    std::vector<Atom> m_fmole1;
    std::vector<Atom> m_xmole2;
  };

  class Molecule
  {
  public:
    explicit Molecule(){}
    void addAtom(std::string& line);

    std::vector<Atom> m_atoms;
    double m_energy; //energe
  };

  class DataBase
  {
  };

  class EnergeForceDatabase
  {
  public:
    explicit EnergeForceDatabase(const char* filename, const char* fragtype);
    ~EnergeForceDatabase();
    Config*& at_all(int dim1, int dim2, int dim3, int dim4);
    Config*& at_wtr(int dim1, int dim2);
    void set_wtr(int dim1, int dim2, Config* cfg);
    void set_all(int dim1, int dim2, int dim3, int dim4, Config* cfg);
    void read_file();

  private:

    // read_mole from the first two line
    void read_mole(int ind, igzstream& datafile);
    double H_energy(Config& config);
    double charge(double q0, double q1, double r);

    Config** m_all_config;
    Config** m_wtr_config;
    std::string m_filename;
    std::string m_fragtype;
    DataStructure* m_pDS;
    int m_natoms1;
    int m_natoms2;

    Molecule m_mole1;
    Molecule m_mole2;

  };
}

#endif
